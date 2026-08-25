//! Calibration tables for a `.d`, optionally via Bruker's `timsdata`.
//!
//! Nothing Bruker-derived is compiled into or shipped with this crate. The
//! library is located at *runtime* from a path the caller supplies, loaded with
//! `libloading`, and used for exactly one purpose: building the two calibration
//! lookup tables for a `.d`.
//!
//! That is the whole of what the SDK contributes to centroiding. The
//! centroiding algorithm itself works in tof/scan index space and never calls
//! into it. So the SDK is consulted once per file, up front, to answer "what
//! m/z is tof index N" and "what 1/K0 is scan number N" across the full index
//! range; afterwards conversion is an array index and the library is closed.
//!
//! When no path is supplied, callers fall back to timsrust's own converters,
//! which derive an approximate calibration from the acquisition range in
//! `analysis.tdf`. Those are *not* equivalent to the SDK's values.

use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::path::Path;

use libloading::{Library, Symbol};

/// Signature of `tims_open_v2`.
type TimsOpenV2 = unsafe extern "C" fn(*const c_char, u32, u32) -> u64;
/// Signature of `tims_close`.
type TimsClose = unsafe extern "C" fn(u64);
/// Shared signature of the four index/value conversion entry points.
type TimsConvert = unsafe extern "C" fn(u64, i64, *const f64, *mut f64, u32) -> u32;
/// Signature of `tims_get_last_error_string`.
type TimsLastError = unsafe extern "C" fn(*mut c_char, u32) -> u32;

/// Bruker's own default: use the acquisition-time calibration rather than a
/// later DataAnalysis recalibration. Matches what timsrust-sdk requests, so
/// output stays comparable with the stock centroider.
const USE_RECALIBRATED_STATE: u32 = 0;
/// No pressure compensation, again matching timsrust-sdk's behaviour.
const PRESSURE_COMPENSATION_NONE: u32 = 0;

#[derive(Debug)]
pub enum SdkError {
    Load(String),
    Symbol(String),
    Open(String),
    Convert(String),
}

impl std::fmt::Display for SdkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SdkError::Load(m) => write!(f, "could not load Bruker timsdata library: {m}"),
            SdkError::Symbol(m) => write!(f, "missing symbol in Bruker timsdata library: {m}"),
            SdkError::Open(m) => write!(f, "Bruker SDK could not open the .d folder: {m}"),
            SdkError::Convert(m) => write!(f, "Bruker SDK conversion failed: {m}"),
        }
    }
}

impl std::error::Error for SdkError {}

/// The two dense calibration tables for one `.d`, as produced by the SDK.
///
/// Indexed directly: `mz[tof_index]` and `inv_ion_mobility[scan_number]`.
pub struct CalibrationTables {
    pub mz: Vec<f64>,
    pub inv_ion_mobility: Vec<f64>,
}

/// An opened Bruker library plus a handle to one `.d`.
///
/// The `Library` must outlive every symbol resolved from it, so it is kept
/// alongside the handle and both are released together in `Drop`.
pub struct TimsDataLib {
    lib: Library,
    handle: u64,
}

impl TimsDataLib {
    /// Load `timsdata` from `sdk_path` and open `d_path`.
    ///
    /// `sdk_path` may be the library file itself (`libtimsdata.so`,
    /// `timsdata.dll`) or a directory containing it.
    pub fn open(sdk_path: &str, d_path: &str) -> Result<Self, SdkError> {
        let lib_file = resolve_library_path(sdk_path)?;

        // SAFETY: loading an arbitrary shared library runs its initialisers.
        // The caller is explicitly pointing at a Bruker SDK install.
        let lib = unsafe { Library::new(&lib_file) }
            .map_err(|e| SdkError::Load(format!("{} ({e})", lib_file.display())))?;

        let handle = {
            let open: Symbol<TimsOpenV2> = unsafe { lib.get(b"tims_open_v2\0") }
                .map_err(|e| SdkError::Symbol(format!("tims_open_v2 ({e})")))?;
            let c_path = CString::new(d_path)
                .map_err(|e| SdkError::Open(format!("path contains a NUL byte ({e})")))?;
            unsafe { open(c_path.as_ptr(), USE_RECALIBRATED_STATE, PRESSURE_COMPENSATION_NONE) }
        };

        if handle == 0 {
            let msg = last_error(&lib).unwrap_or_else(|| "tims_open_v2 returned 0".to_string());
            return Err(SdkError::Open(msg));
        }

        Ok(Self { lib, handle })
    }

    /// Build both calibration tables in two batched SDK calls.
    ///
    /// `tof_max_index` is `DigitizerNumSamples` and `scan_max_index` the largest
    /// `NumScans` across frames; both tables are inclusive of their maximum, so
    /// any valid index from the raw data is in range.
    ///
    /// Frame 1 is used for both conversions, matching timsrust-sdk. That is
    /// correct when a run carries a single calibration, which is the usual case
    /// -- `Frames.MzCalibration` and `Frames.TimsCalibration` are foreign keys,
    /// so a file with more than one row in either table would need per-frame
    /// tables instead.
    pub fn calibration_tables(
        &self,
        tof_max_index: u32,
        scan_max_index: u32,
    ) -> Result<CalibrationTables, SdkError> {
        let index_to_mz: Symbol<TimsConvert> = unsafe { self.lib.get(b"tims_index_to_mz\0") }
            .map_err(|e| SdkError::Symbol(format!("tims_index_to_mz ({e})")))?;
        let scan_to_im: Symbol<TimsConvert> =
            unsafe { self.lib.get(b"tims_scannum_to_oneoverk0\0") }
                .map_err(|e| SdkError::Symbol(format!("tims_scannum_to_oneoverk0 ({e})")))?;

        let mz = self.batch(&index_to_mz, tof_max_index, "tims_index_to_mz")?;
        let inv_ion_mobility = self.batch(&scan_to_im, scan_max_index, "tims_scannum_to_oneoverk0")?;

        Ok(CalibrationTables { mz, inv_ion_mobility })
    }

    /// Evaluate one conversion across `0..=max_index` in a single call.
    fn batch(
        &self,
        func: &Symbol<TimsConvert>,
        max_index: u32,
        name: &str,
    ) -> Result<Vec<f64>, SdkError> {
        let n = max_index as usize + 1;
        let input: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let mut output = vec![0.0f64; n];

        // The C API takes the count as u32, so a single call covers any
        // plausible index range (tof tops out around 10^6).
        let rc = unsafe {
            func(
                self.handle,
                1,
                input.as_ptr(),
                output.as_mut_ptr(),
                n as u32,
            )
        };
        if rc == 0 {
            let msg = last_error(&self.lib).unwrap_or_else(|| format!("{name} returned 0"));
            return Err(SdkError::Convert(msg));
        }
        Ok(output)
    }
}

impl Drop for TimsDataLib {
    fn drop(&mut self) {
        if self.handle != 0 {
            if let Ok(close) = unsafe { self.lib.get::<TimsClose>(b"tims_close\0") } {
                unsafe { close(self.handle) };
            }
        }
    }
}

/// Accept either the library file or a directory containing it.
fn resolve_library_path(sdk_path: &str) -> Result<std::path::PathBuf, SdkError> {
    let p = Path::new(sdk_path);
    if p.is_file() {
        return Ok(p.to_path_buf());
    }
    if p.is_dir() {
        // Platform-specific names, tried in order; the caller's platform
        // determines which one exists.
        for name in ["libtimsdata.so", "timsdata.dll", "libtimsdata.dylib"] {
            let candidate = p.join(name);
            if candidate.is_file() {
                return Ok(candidate);
            }
        }
        return Err(SdkError::Load(format!(
            "no libtimsdata.so / timsdata.dll found in directory {sdk_path}"
        )));
    }
    Err(SdkError::Load(format!("path does not exist: {sdk_path}")))
}

/// Pull the SDK's thread-local error string, if it has one.
fn last_error(lib: &Library) -> Option<String> {
    let f: Symbol<TimsLastError> = unsafe { lib.get(b"tims_get_last_error_string\0") }.ok()?;
    let mut buf = vec![0i8; 2048];
    let len = unsafe { f(buf.as_mut_ptr() as *mut c_char, buf.len() as u32) };
    if len == 0 {
        return None;
    }
    let s = unsafe { CStr::from_ptr(buf.as_ptr() as *const c_char) };
    Some(s.to_string_lossy().into_owned())
}
