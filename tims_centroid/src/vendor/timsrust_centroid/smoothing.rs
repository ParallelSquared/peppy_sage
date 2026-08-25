use rustc_hash::FxHashMap;

#[derive(Debug, Clone)]
pub(crate) struct ScanSmoother {
    scan_kernel: Vec<u64>,
    min_count: usize,
}

impl ScanSmoother {
    pub(crate) fn new(scan_kernel: Vec<u64>, min_count: usize) -> Self {
        Self {
            scan_kernel,
            min_count,
        }
    }

    pub(crate) fn kernel(&self) -> &[u64] {
        &self.scan_kernel
    }

    pub(crate) fn len(&self) -> usize {
        self.scan_kernel.len()
    }

    pub(crate) fn min_count(&self) -> usize {
        self.min_count
    }

    /// Smooth `options` along the scan axis, reusing `scratch`.
    ///
    /// MODIFIED FROM UPSTREAM. The original allocated three containers on every
    /// call -- the accumulator map, the sorted scan list, and the usable-flag
    /// vector -- and this runs once per tof index per frame, so on a
    /// 67,778-frame acquisition that was the single largest source of
    /// allocator traffic. They are now caller-owned buffers, cleared and
    /// refilled in place so their capacity survives across calls.
    ///
    /// The returned map is still freshly allocated: `RollingSparseBuffer::
    /// rollover` takes ownership of it, so it cannot come from scratch.
    ///
    /// Behaviour is unchanged.
    pub(crate) fn smooth(
        &self,
        options: &FxHashMap<u32, (u64, u32)>,
        scratch: &mut SmoothScratch,
    ) -> FxHashMap<u32, u64> {
        // Destructured so the three buffers can be borrowed independently.
        let SmoothScratch { opt, scans, usable } = scratch;
        opt.clear();
        scans.clear();
        usable.clear();

        scans.extend(options.keys().copied());
        scans.sort_unstable();
        usable.resize(scans.len(), false);
        Self::set_usable_candidates(
            self.scan_kernel.len() as u32,
            options,
            self.min_count,
            scans,
            usable,
        );

        for (i, &scan) in scans.iter().enumerate() {
            if !usable[i] {
                continue;
            }
            let &(apex_intensity, count) = options.get(&scan).unwrap();
            for (scan_offset, &value) in self.scan_kernel.iter().enumerate() {
                let key = scan + scan_offset as u32;
                let entry = opt.entry(key).or_insert((0, 0));
                entry.0 += apex_intensity * value;
                entry.1 += count;
            }
        }
        opt.iter()
            .filter_map(|(&scan, &(apex_intensity, count))| {
                if count >= self.min_count as u32 {
                    Some((scan, apex_intensity))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Mark which scans sit in a window carrying enough ions to be kept.
    ///
    /// MODIFIED FROM UPSTREAM: writes into a caller-supplied buffer, which the
    /// caller has already sized and zeroed, rather than returning a fresh
    /// `Vec<bool>`.
    fn set_usable_candidates(
        kernel_len: u32,
        options: &FxHashMap<u32, (u64, u32)>,
        min_count: usize,
        scans: &[u32],
        result: &mut [bool],
    ) {
        let mut left = 0;
        let mut right = 0;
        while left < scans.len() {
            while right < scans.len() && scans[right] < scans[left] + kernel_len
            {
                right += 1;
            }
            let window_count: u32 = scans[left..right]
                .iter()
                .map(|scan| options.get(scan).unwrap().1)
                .sum();
            if window_count >= min_count as u32 {
                for i in result.iter_mut().take(right).skip(left) {
                    *i = true;
                }
            }
            left += 1;
        }
    }
}

/// Reusable working buffers for [`ScanSmoother::smooth`].
///
/// Held by the caller across calls so the allocations amortise. `clear` keeps
/// capacity, so after the first few tof indices the smoother stops allocating
/// almost entirely.
#[derive(Debug, Clone, Default)]
pub(crate) struct SmoothScratch {
    opt: FxHashMap<u32, (u64, u32)>,
    scans: Vec<u32>,
    usable: Vec<bool>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smooth_basic() {
        let smoother = ScanSmoother::new(vec![1, 2, 1], 1);
        let mut options = FxHashMap::default();
        options.insert(1, (10, 1));
        options.insert(2, (20, 1));
        options.insert(4, (40, 1));
        let result =
            smoother.smooth(&options, &mut SmoothScratch::default());
        let mut expected = FxHashMap::default();
        expected.insert(1, 10);
        expected.insert(2, 40);
        expected.insert(3, 50);
        expected.insert(4, 60);
        expected.insert(5, 80);
        expected.insert(6, 40);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_smooth_min_count() {
        let smoother = ScanSmoother::new(vec![1, 1], 3);
        let mut options = FxHashMap::default();
        options.insert(1, (5, 1));
        options.insert(2, (10, 1));
        options.insert(3, (15, 1));
        // Each scan has count 1, so any window of size 2 has at most 2 counts
        // min_count is 3, so result should be empty
        let result =
            smoother.smooth(&options, &mut SmoothScratch::default());
        assert!(result.is_empty());
    }

    #[test]
    fn test_set_usable_candidates() {
        let mut options = FxHashMap::default();
        options.insert(1, (10, 1));
        options.insert(2, (20, 2));
        options.insert(4, (40, 1));
        let scans = vec![1u32, 2, 4];
        let mut usable = vec![false; scans.len()];
        ScanSmoother::set_usable_candidates(
            2,
            &options,
            3,
            &scans,
            &mut usable,
        );
        assert_eq!(usable, vec![true, true, false]);
    }
}
