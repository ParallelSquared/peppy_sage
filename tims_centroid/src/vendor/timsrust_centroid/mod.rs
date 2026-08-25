//! 2D centroiding for Bruker timsTOF data.
//!
//! # Provenance and licence
//!
//! Vendored from the `timsrust-centroid` crate at version 0.6.4
//! (<https://github.com/mannlabs/timsrust>), which is licensed under the
//! **Apache License, Version 2.0** -- see `LICENSE-APACHE` beside this module.
//! The surrounding `peppy_sage` crate is MIT; that difference is why this code
//! keeps its own licence file rather than inheriting the crate's.
//!
//! It is vendored rather than depended upon because the upstream algorithm
//! allocates per frame in a way that dominates runtime (see below), and the fix
//! has to ship inside a published `peppy_sage` wheel. A path or git dependency
//! could not be published.
//!
//! `spectrum_reader` is deliberately not vendored: the parquet path only needs
//! `PeakReader`.
//!
//! # Modifications from upstream
//!
//! Profiling a full 67,778-frame run put ~55-60% of active CPU in `malloc`,
//! `free`, `memmove` and hashmap rehashing rather than in the centroiding
//! arithmetic, because `ScanSmoother::smooth` allocated three fresh containers
//! on every call and is invoked once per tof index per frame. The changes here
//! reuse that scratch space across calls instead. See `smoothing.rs`.

mod buffer;
mod centroider;
mod error;
mod peakbuffer;
mod peaks;
mod smoothing;

pub use error::{TimsError, TimsResult};
pub use peaks::{Peak, PeakReader};
