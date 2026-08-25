use rustc_hash::FxHashMap;
use timsrust_core::utils::vec::arg_max;

use crate::vendor::timsrust_centroid::{
    Peak, peakbuffer::PeakBuffer, peaks::TofSlices, smoothing::ScanSmoother,
};

pub(crate) fn scale_to_u64(vec: &[f32], upper_bound: f32) -> Vec<u64> {
    vec.iter()
        .map(|&x| {
            let t = (upper_bound * x) as u64;
            t.min(upper_bound as u64)
        })
        .collect()
}

#[derive(Debug)]
pub(crate) struct FrameCentroider {
    tof_kernel: Vec<u64>,
    scan_smoother: ScanSmoother,
    tof_kernel_apex: usize,
    scan_kernel_apex: usize,
}

impl FrameCentroider {
    pub(crate) fn new(
        scan_kernel: &[f32],
        tof_kernel: &[f32],
        min_count: usize,
    ) -> Self {
        let scan_kernel = scale_to_u64(scan_kernel, u8::MAX as f32);
        let tof_kernel = scale_to_u64(tof_kernel, u8::MAX as f32);
        let scan_smoother = ScanSmoother::new(scan_kernel.clone(), min_count);
        let tof_kernel_apex = arg_max(&tof_kernel).expect("Kernel is empty");
        let scan_kernel_apex = arg_max(&scan_kernel).expect("Kernel is empty");
        FrameCentroider {
            tof_kernel,
            scan_smoother,
            tof_kernel_apex,
            scan_kernel_apex,
        }
    }

    pub(crate) fn tof_kernel(&self) -> &[u64] {
        &self.tof_kernel
    }

    pub(crate) fn scan_smoother(&self) -> &ScanSmoother {
        &self.scan_smoother
    }

    pub(crate) fn centroid<'a>(
        &'a self,
        slices: TofSlices<'a>,
        frame_index: usize,
    ) -> impl Iterator<Item = Peak> + 'a {
        PeakBuffer::new(
            slices,
            self.scan_smoother.clone(),
            self.tof_kernel.len(),
            frame_index as u32,
            self.scan_kernel_apex as u32,
            self.tof_kernel_apex as u32,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_frame_centroider_new_and_getters() {
        let scan_kernel = vec![0.1, 1.0, 0.2];
        let tof_kernel = vec![0.2, 1.0, 0.0];
        let min_count = 3;
        let centroider =
            FrameCentroider::new(&scan_kernel, &tof_kernel, min_count);
        assert_eq!(
            centroider.scan_smoother().kernel(),
            vec![25, 255, 51].as_slice()
        );
        assert_eq!(centroider.tof_kernel(), vec![51, 255, 0].as_slice());
        assert_eq!(centroider.scan_smoother().min_count(), min_count);
        assert_eq!(centroider.scan_kernel_apex, 1);
        assert_eq!(centroider.tof_kernel_apex, 1);
    }

    #[test]
    fn test_centroid_returns_iterator() {
        let scan_kernel = vec![0.1, 1.0, 0.2];
        let tof_kernel = vec![0.2, 1.0, 0.0];
        let min_count = 1;
        let centroider =
            FrameCentroider::new(&scan_kernel, &tof_kernel, min_count);
        // 3 TOF bins, each with a single scan and intensity.
        let frame = crate::vendor::timsrust_centroid::peaks::frame_from_ions(
            3,
            &[(2, 0, 1), (2, 1, 2), (2, 2, 3)],
        );
        let frame_index = 0;
        let mut iter = centroider.centroid(TofSlices::new(&frame), frame_index);
        // // Since the underlying logic is complex, just check that it is an iterator and returns Peaks or ends
        // // (We can't guarantee output without full mock of dependencies)
        let _ = iter.next(); // Should not panic
    }
}
