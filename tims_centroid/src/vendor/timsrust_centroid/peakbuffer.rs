use rustc_hash::FxHashMap;

use crate::vendor::timsrust_centroid::{
    Peak, buffer::RollingSparseBuffer, peaks::TofSlices,
    smoothing::{ScanSmoother, SmoothScratch},
};

type Options = FxHashMap<u32, (u64, u32)>;

#[derive(Debug)]
pub(crate) struct PeakBuffer<'a> {
    buffer: RollingSparseBuffer,
    options: Vec<Options>,
    slices: TofSlices<'a>,
    /// Ring of the `tof_kernel_len + 1` tof slices currently in view, indexed
    /// by `tof_index % window.len()`.
    window: Vec<Vec<(u32, u64)>>,
    prev_hits: FxHashMap<u64, Vec<Peak>>,
    new_hits: FxHashMap<u64, Vec<Peak>>,
    peak_queue: Vec<Peak>,
    scan_smoother: ScanSmoother,
    /// Reused across every smooth() call for this frame; see SmoothScratch.
    smooth_scratch: SmoothScratch,
    tof_kernel_len: usize,
    frame_index: u32,
    total_len: usize,
    scan_kernel_apex: u32,
    tof_kernel_apex: u32,
    tof_index: usize,
}

impl<'a> PeakBuffer<'a> {
    pub(crate) fn new(
        mut slices: TofSlices<'a>,
        scan_smoother: ScanSmoother,
        tof_kernel_len: usize,
        frame_index: u32,
        scan_kernel_apex: u32,
        tof_kernel_apex: u32,
    ) -> Self {
        assert!(tof_kernel_len >= 3);
        let max_scan = slices.max_scan()
            + 1
            + scan_kernel_apex
            + scan_smoother.len() as u32;
        let mut buffer = RollingSparseBuffer::new(max_scan);
        // Indices past the last one carrying signal are empty, exactly as the
        // trailing empty maps used to be.
        let total_len = slices.max_tof() + buffer.len();

        // Prime the window with slices 0..=tof_kernel_len; the first rollover
        // reads index tof_kernel_len and index 0 in the same step.
        let mut window: Vec<Vec<(u32, u64)>> =
            vec![Vec::new(); tof_kernel_len + 1];
        for t in 0..=tof_kernel_len {
            let slot = t % window.len();
            slices.fill(t, &mut window[slot]);
        }
        let mut options = Vec::with_capacity(tof_kernel_len);
        let mut first = Options::default();
        Self::add_to_options(&window[0], &mut first);
        options.push(first);
        for t in 1..tof_kernel_len {
            let mut new_options =
                options.last().expect("No last option").clone();
            Self::add_to_options(&window[t % (tof_kernel_len + 1)], &mut new_options);
            options.push(new_options);
        }
        let mut smooth_scratch = SmoothScratch::default();
        options.iter().take(buffer.len()).for_each(|opt| {
            let opt = scan_smoother.smooth(opt, &mut smooth_scratch);
            buffer.rollover(opt);
        });
        let prev_hits: FxHashMap<u64, Vec<Peak>> = FxHashMap::default();
        let new_hits: FxHashMap<u64, Vec<Peak>> = FxHashMap::default();
        let peak_queue: Vec<Peak> = Vec::new();
        Self {
            buffer,
            options,
            slices,
            window,
            prev_hits,
            new_hits,
            peak_queue,
            scan_smoother,
            smooth_scratch,
            tof_kernel_len,
            frame_index,
            total_len,
            scan_kernel_apex,
            tof_kernel_apex,
            tof_index: 0,
        }
    }

    pub(crate) fn rollover(&mut self) {
        self.collect_peaks();
        // Recycle the window entry falling off the front rather than cloning
        // the last one into a fresh allocation; clone_from reuses its table.
        let mut new_options = self.options.remove(0);
        new_options.clone_from(self.options.last().expect("No last option"));
        let n = self.window.len();
        Self::add_to_options(
            &self.window[(self.tof_index + self.tof_kernel_len) % n],
            &mut new_options,
        );
        Self::remove_from_options(
            &self.window[self.tof_index % n],
            &mut new_options,
        );
        self.options.push(new_options);
        let opt = self
            .scan_smoother
            .smooth(&self.options[self.buffer.len() - 1], &mut self.smooth_scratch);
        self.buffer.rollover(opt);
        self.tof_index += 1;
        // The slice that just left the window is reused for the one entering it.
        let n = self.window.len();
        let entering = self.tof_index + self.tof_kernel_len;
        let slot = entering % n;
        let mut buf = std::mem::take(&mut self.window[slot]);
        self.slices.fill(entering, &mut buf);
        self.window[slot] = buf;
    }

    fn remove_from_options(
        tofs: &[(u32, u64)],
        new_options: &mut Options,
    ) {
        for &(scan, apex_intensity) in tofs.iter() {
            let scan_index = scan;
            let entry = new_options.entry(scan_index).or_insert((0, 0));
            entry.0 -= apex_intensity;
            entry.1 -= 1;
            if entry.1 == 0 {
                new_options.remove(&scan_index);
            }
        }
    }

    fn add_to_options(tofs: &[(u32, u64)], new_options: &mut Options) {
        for &(scan, apex_intensity) in tofs.iter() {
            let scan_index = scan;
            let entry = new_options.entry(scan_index).or_insert((0, 0));
            entry.0 += apex_intensity;
            entry.1 += 1;
        }
    }

    fn collect_peaks(&mut self) {
        for (scan_index, apex_intensity, unique_peak) in
            self.buffer.iter_peaks(self.scan_smoother.len() / 2)
        {
            let peak = Peak {
                frame: self.frame_index,
                scan: scan_index - self.scan_kernel_apex,
                tof: 1 + self.tof_index as u32 - self.tof_kernel_apex,
                apex_intensity,
            };
            if unique_peak {
                self.peak_queue.push(peak);
            } else {
                self.new_hits.entry(apex_intensity).or_default().push(peak);
            }
        }
        for (apex_intensity, ambiguous_peaks) in self.prev_hits.drain() {
            match self.new_hits.get_mut(&apex_intensity) {
                Some(existing_peaks) => {
                    existing_peaks.extend(ambiguous_peaks);
                },
                None => {
                    let peak = Self::pick(ambiguous_peaks);
                    self.peak_queue.push(peak);
                },
            }
        }
        std::mem::swap(&mut self.prev_hits, &mut self.new_hits);
    }

    fn pick(ambiguous_peaks: Vec<Peak>) -> Peak {
        // TODO
        let len = ambiguous_peaks.len() as u32;
        let mut sum = Peak::default();
        for peak in ambiguous_peaks.iter() {
            sum.frame += peak.frame;
            sum.scan += peak.scan;
            sum.tof += peak.tof;
            sum.apex_intensity += peak.apex_intensity;
        }
        Peak {
            frame: sum.frame / len,
            scan: sum.scan / len,
            tof: sum.tof / len,
            apex_intensity: sum.apex_intensity / len as u64,
        }
    }
}

impl Iterator for PeakBuffer<'_> {
    type Item = Peak;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(peak) = self.peak_queue.pop() {
                return Some(peak);
            }
            if self.tof_index > (self.total_len - self.tof_kernel_len - 1) {
                return None;
            }
            self.rollover();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vendor::timsrust_centroid::peaks::frame_from_ions;

    #[test]
    fn test_peakbuffer_new_and_iter() {
        // (scan, tof, intensity); same signal the old Vec<HashMap> fixture held.
        let frame = frame_from_ions(
            4,
            &[
                (1, 0, 100), (2, 0, 200),
                (1, 1, 150), (3, 1, 250),
                (2, 2, 120),
                (1, 3, 180),
            ],
        );
        let tofs = TofSlices::new(&frame);
        let scan_kernel = vec![1, 2, 1];
        let scan_smoother = ScanSmoother::new(scan_kernel, 1);
        let tof_len = 3;
        let frame_index = 0;
        let scan_kernel_apex = 1;
        let tof_kernel_apex = 1;
        let mut pb = PeakBuffer::new(
            tofs,
            scan_smoother,
            tof_len,
            frame_index,
            scan_kernel_apex,
            tof_kernel_apex,
        );
        let _ = pb.next();
    }
}
