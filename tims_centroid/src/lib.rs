//! 2D centroiding of timsTOF `.d` data to a `peaks.parquet`.
//!
//! This mirrors the `run_parquet` path of `timsrust-centroid-cli` -- same
//! algorithm, same output schema -- with one difference: calibration is
//! supplied at runtime rather than linked at build time.
//!
//! The stock CLI reaches Bruker's calibration through a Cargo feature that
//! links `timsdata` at build time, which means shipping Bruker's library
//! alongside the binary. Here the centroiding is done by `timsrust-centroid` in
//! tof/scan index space (it needs no calibration at all), and m/z and 1/K0 are
//! filled in at the end from either:
//!
//! - lookup tables built by the Bruker SDK, loaded at runtime from a
//!   caller-supplied path (see [`crate::calibration`]), or
//! - timsrust's own converters, which approximate the calibration from the
//!   acquisition range in `analysis.tdf`, when no SDK path is given.
//!
//! Retention time never comes from the SDK -- it is read from the `Frames`
//! table either way.

pub mod calibration;
mod vendor;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use timsrust::{ImConverter, MzConverter, RtConverter};
use crate::vendor::timsrust_centroid::{PeakReader, TimsError, TimsResult};
use timsrust_core::io::formats::parquet::ParquetWriter;
use timsrust_core::utils::thread::Synced;
use timsrust_core::{ConvertibleTo, FrameIndex, ScanIndex, TofIndex};
use timsrust_tdf::TdfFrameReader;

use crate::calibration::{CalibrationTables, TimsDataLib};

/// One centroided peak. Field order and types match the stock centroider's
/// `peaks.parquet` so downstream readers need no changes.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FullPeak {
    pub frame: u32,
    pub scan: u32,
    pub tof: u32,
    pub apex_intensity: u64,
    pub rt: f64,
    pub im: f64,
    pub mz: f64,
    pub isolation_window_lower: Option<f64>,
    pub isolation_window_upper: Option<f64>,
}

/// How tof/scan indices become m/z and 1/K0.
///
/// Both variants are pure lookups by the time centroiding output is written;
/// the SDK, when used, has already been consulted and closed.
enum Calibration {
    Sdk(CalibrationTables),
    Estimated {
        mz: MzConverter,
        im: ImConverter,
    },
}

impl Calibration {
    fn mz(&self, tof: u32) -> f64 {
        match self {
            // Indices come from the instrument and are always within the
            // digitizer range; clamp so a malformed frame cannot panic.
            Calibration::Sdk(t) => t.mz[(tof as usize).min(t.mz.len() - 1)],
            Calibration::Estimated { mz, .. } => f64::from(
                TofIndex::try_from(tof)
                    .expect("tof index out of range")
                    .convert(mz),
            ),
        }
    }

    fn im(&self, scan: u32) -> f64 {
        match self {
            Calibration::Sdk(t) => {
                t.inv_ion_mobility[(scan as usize).min(t.inv_ion_mobility.len() - 1)]
            }
            Calibration::Estimated { im, .. } => f64::from(
                ScanIndex::try_from(scan)
                    .expect("scan index out of range")
                    .convert(im),
            ),
        }
    }
}

/// Aim for roughly this many progress updates over a whole run.
const PROGRESS_UPDATES: usize = 100;
/// Never report more often than this, however short the run.
const PROGRESS_INTERVAL_MIN: usize = 1;
/// Never report less often than this, however long the run.
///
/// This also bounds how quickly a cancellation is noticed, since the callback
/// is the only place the run checks whether it should stop -- at ~120 frames/s
/// this keeps Ctrl-C responding in about a second, at a cost of a few hundred
/// GIL acquisitions across a full run.
const PROGRESS_INTERVAL_MAX: usize = 128;

/// Frames between progress callbacks, scaled to the size of the run.
///
/// Each callback reacquires the GIL from a rayon worker, so a fixed small
/// interval would serialise the pool on a long run. A fixed large one would
/// leave a short run with a single update at the end. Scaling keeps the bar
/// informative either way.
fn progress_interval(frame_count: usize) -> usize {
    (frame_count / PROGRESS_UPDATES)
        .clamp(PROGRESS_INTERVAL_MIN, PROGRESS_INTERVAL_MAX)
}

/// Centroid `in_path` (a `.d` folder) and write `out_path` as parquet.
///
/// `sdk_path` points at Bruker's `timsdata` library, or a directory containing
/// it. When `None`, the approximate calibration is used and the m/z and 1/K0
/// columns will *not* match an SDK-calibrated run.
///
/// `progress` is called with `(frames_done, frames_total)` every
/// `progress_interval(frame_count)` frames, and once more on completion. It is called from
/// rayon worker threads, so it must be cheap and must not assume an order.
///
/// Returning `false` from `progress` cancels the run: remaining frames are
/// skipped and the function returns early. The caller owns the partially
/// written output. This is how an interrupt reaches a loop that would otherwise
/// run to completion with the GIL released.
///
/// Returns the number of peaks written.
pub fn centroid_to_parquet(
    in_path: &str,
    out_path: &str,
    min_ion_count_ms1: f64,
    min_ion_count_ms2: f64,
    sdk_path: Option<&str>,
    progress: Option<&(dyn Fn(usize, usize) -> bool + Sync)>,
) -> TimsResult<usize> {
    let frame_reader = TdfFrameReader::new(in_path)
        .map(|r| r.into_inner())
        .map_err(|e| TimsError::new(e.to_string()))?;
    let peak_reader = PeakReader::new(frame_reader, min_ion_count_ms1, min_ion_count_ms2)?;

    let rt_converter = RtConverter::new(in_path)
        .ok_or_else(|| TimsError::new("could not build retention time converter".to_string()))?;

    let calibration = match sdk_path {
        Some(path) => {
            let (tof_max_index, scan_max_index) = index_bounds(in_path)?;
            let sdk = TimsDataLib::open(path, in_path).map_err(|e| TimsError::new(e.to_string()))?;
            let tables = sdk
                .calibration_tables(tof_max_index, scan_max_index)
                .map_err(|e| TimsError::new(e.to_string()))?;
            // Tables are built; the library handle closes here.
            Calibration::Sdk(tables)
        },
        None => Calibration::Estimated {
            mz: MzConverter::new(in_path)
                .ok_or_else(|| TimsError::new("could not build m/z converter".to_string()))?,
            im: ImConverter::new(in_path)
                .ok_or_else(|| TimsError::new("could not build mobility converter".to_string()))?,
        },
    };

    let frame_count = peak_reader.frame_count();
    let mut writer = ParquetWriter::<FullPeak>::new(out_path)
        .map_err(|e| TimsError::new(e.to_string()))?;
    writer.set_max_batch_count(frame_count);
    let synced_writer = Synced::from(writer);

    // Frames complete out of order, so progress is a count of finished frames
    // rather than the highest index reached.
    let frames_done = std::sync::atomic::AtomicUsize::new(0);
    let interval = progress_interval(frame_count);

    // rayon cannot break out of a parallel iterator, so cancellation makes the
    // remaining frames no-ops. Skipping is cheap next to centroiding one.
    let cancelled = std::sync::atomic::AtomicBool::new(false);

    let peak_count = (0..frame_count)
        .into_par_iter()
        .map(|index| {
            if cancelled.load(std::sync::atomic::Ordering::Relaxed) {
                return 0;
            }

            let result = centroid_one_frame(
                index,
                &peak_reader,
                &calibration,
                &rt_converter,
                &synced_writer,
            );

            if let Some(cb) = progress {
                let done = frames_done
                    .fetch_add(1, std::sync::atomic::Ordering::Relaxed)
                    + 1;
                if done % interval == 0 || done == frame_count {
                    if !cb(done, frame_count) {
                        cancelled.store(true, std::sync::atomic::Ordering::Relaxed);
                    }
                }
            }
            result
        })
        .sum::<usize>();

    Ok(peak_count)
}

/// Centroid one frame and write its peaks. Returns the number written.
fn centroid_one_frame(
    index: usize,
    peak_reader: &PeakReader<
        timsrust_tdf::TdfIonReader,
        timsrust_tdf::FrameInfoReader,
    >,
    calibration: &Calibration,
    rt_converter: &RtConverter,
    synced_writer: &Synced<ParquetWriter<FullPeak>>,
) -> usize {
    {
            let Ok(peaks) = peak_reader.get_peaks_from_frame(index) else {
                return 0;
            };
            // Only the frame's quadrupole settings are needed here, so read the
            // info without decoding the frame's ion data again.
            let Ok(frame) = peak_reader
                .frame_reader()
                .get_partial_frame_without_ions(index)
            else {
                return 0;
            };
            let quad_info = frame.info().quadrupole_settings();

            // Retention time is a property of the frame, not the peak.
            let rt = peaks
                .first()
                .map(|p| {
                    f64::from(
                        FrameIndex::try_from(p.frame)
                            .expect("frame index out of range")
                            .convert(&rt_converter),
                    )
                })
                .unwrap_or(0.0);

            let peaks: Vec<FullPeak> = peaks
                .into_iter()
                .map(|p| {
                    let window = quad_info.get_isolation_window(p.scan as usize);
                    FullPeak {
                        frame: p.frame,
                        scan: p.scan,
                        tof: p.tof,
                        apex_intensity: p.apex_intensity,
                        rt,
                        im: calibration.im(p.scan),
                        mz: calibration.mz(p.tof),
                        isolation_window_lower: window
                            .as_ref()
                            .map(|w| f64::from(w.lower())),
                        isolation_window_upper: window
                            .as_ref()
                            .map(|w| f64::from(w.upper())),
                    }
                })
                .collect();

            let written = peaks.len();
            let _ = synced_writer
                .with_lock(|w| w.write_batch(peaks))
                .expect("parquet writer mutex poisoned");
            written
    }
}

/// Largest valid tof and scan indices, for sizing the SDK lookup tables.
///
/// Read straight from `analysis.tdf`: `DigitizerNumSamples` bounds the tof axis
/// and the largest per-frame `NumScans` bounds the mobility axis.
fn index_bounds(d_path: &str) -> TimsResult<(u32, u32)> {
    let tdf = std::path::Path::new(d_path).join("analysis.tdf");
    let conn = rusqlite::Connection::open(&tdf)
        .map_err(|e| TimsError::new(format!("could not open {}: {e}", tdf.display())))?;

    let digitizer: String = conn
        .query_row(
            "SELECT Value FROM GlobalMetadata WHERE Key = 'DigitizerNumSamples'",
            [],
            |row| row.get(0),
        )
        .map_err(|e| TimsError::new(format!("DigitizerNumSamples not readable: {e}")))?;
    let tof_max_index: u32 = digitizer
        .trim()
        .parse()
        .map_err(|e| TimsError::new(format!("DigitizerNumSamples is not an integer: {e}")))?;

    let scan_max_index: u32 = conn
        .query_row("SELECT MAX(NumScans) FROM Frames", [], |row| row.get(0))
        .map_err(|e| TimsError::new(format!("MAX(NumScans) not readable: {e}")))?;

    Ok((tof_max_index, scan_max_index))
}
