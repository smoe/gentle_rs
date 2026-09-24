//! GC-content computations and display helpers.

use serde::{Deserialize, Serialize};

pub const DEFAULT_SECTION_SIZE_BP: usize = 100;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GcRegion {
    from: usize,
    to: usize,
    gc: f32,
}

impl GcRegion {
    #[inline(always)]
    pub fn from(&self) -> usize {
        self.from
    }

    #[inline(always)]
    pub fn to(&self) -> usize {
        self.to
    }

    #[inline(always)]
    pub fn gc(&self) -> f32 {
        self.gc
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct GcContents {
    regions: Vec<GcRegion>,
}

impl GcContents {
    pub fn new_from_sequence(sequence: &[u8]) -> Self {
        Self::new_from_sequence_with_bin_size(sequence, DEFAULT_SECTION_SIZE_BP)
    }

    pub fn new_from_sequence_with_bin_size(sequence: &[u8], bin_size_bp: usize) -> Self {
        let mut ret = Self::default();
        let mut pos = 0;
        let section_size = Self::get_section_size(sequence, bin_size_bp);
        while pos < sequence.len() {
            let to = sequence.len().min(pos + section_size);
            let gc = Self::calculate_gc(&sequence[pos..to]);
            ret.regions.push(GcRegion { from: pos, to, gc });
            pos += section_size;
        }
        ret
    }

    #[inline(always)]
    pub fn regions(&self) -> &[GcRegion] {
        &self.regions
    }

    /// Count bins intersecting a 0-based, half-open viewport without reading bases.
    ///
    /// `None` counts the full sequence. Zero-sized bins are clamped to one base,
    /// as in GC calculation; empty or reversed viewports contain no bins.
    pub fn region_count_for_viewport(
        sequence_length: usize,
        bin_size_bp: usize,
        viewport: Option<(usize, usize)>,
    ) -> usize {
        let (start, end) = viewport.unwrap_or((0, sequence_length));
        let end = end.min(sequence_length);
        if start >= end {
            return 0;
        }
        let section_size = bin_size_bp.max(1);
        end.div_ceil(section_size) - start / section_size
    }

    #[inline(always)]
    fn get_section_size(sequence: &[u8], bin_size_bp: usize) -> usize {
        sequence.len().min(bin_size_bp.max(1))
    }

    #[inline(always)]
    fn calculate_gc(sequence: &[u8]) -> f32 {
        let gc = sequence
            .iter()
            .map(|c| c.to_ascii_uppercase())
            .filter(|&c| c == b'G' || c == b'C')
            .count() as f32;
        gc / sequence.len() as f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_contents() {
        let sequence = b"AAAGGGTTTCCC";
        let gc_contents = GcContents::new_from_sequence(sequence);
        assert_eq!(gc_contents.regions.len(), 1);
        assert_eq!(
            gc_contents.regions[0],
            GcRegion {
                from: 0,
                to: 12,
                gc: 0.5
            }
        );
    }

    #[test]
    fn test_gc_contents_custom_bin_size() {
        let sequence = b"AAAAGGGGTTTTCCCC";
        let gc_contents = GcContents::new_from_sequence_with_bin_size(sequence, 4);
        assert_eq!(gc_contents.regions.len(), 4);
        assert_eq!(gc_contents.regions[0].gc, 0.0);
        assert_eq!(gc_contents.regions[1].gc, 1.0);
        assert_eq!(gc_contents.regions[2].gc, 0.0);
        assert_eq!(gc_contents.regions[3].gc, 1.0);
    }

    #[test]
    fn test_gc_contents_zero_bin_size_is_clamped() {
        let sequence = b"ATGC";
        let gc_contents = GcContents::new_from_sequence_with_bin_size(sequence, 0);
        assert_eq!(gc_contents.regions.len(), 4);
    }

    #[test]
    fn region_count_matches_materialized_bins_and_half_open_overlap() {
        // Exhaustive small geometry, including partial final bins and outside views.
        for length in 0..=24 {
            let sequence = vec![b'N'; length];
            for bin_size in 0..=length + 2 {
                let gc = GcContents::new_from_sequence_with_bin_size(&sequence, bin_size);
                assert_eq!(
                    GcContents::region_count_for_viewport(length, bin_size, None),
                    gc.regions().len()
                );
                for start in 0..=length + 2 {
                    for end in start + 1..=length + 3 {
                        let expected = gc
                            .regions()
                            .iter()
                            .filter(|region| region.from() < end && region.to() > start)
                            .count();
                        assert_eq!(
                            GcContents::region_count_for_viewport(
                                length,
                                bin_size,
                                Some((start, end))
                            ),
                            expected,
                            "length={length}, bin={bin_size}, viewport={start}..{end}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn region_count_handles_empty_ranges_and_integer_limits() {
        for viewport in [None, Some((0, 0)), Some((0, usize::MAX))] {
            assert_eq!(GcContents::region_count_for_viewport(0, 0, viewport), 0);
        }
        for viewport in [(3, 3), (10, 3), (usize::MAX, usize::MAX)] {
            assert_eq!(
                GcContents::region_count_for_viewport(11, 4, Some(viewport)),
                0
            );
        }
        let length = usize::MAX;
        assert_eq!(
            GcContents::region_count_for_viewport(length, 0, None),
            length
        );
        assert_eq!(
            GcContents::region_count_for_viewport(length, length, None),
            1
        );
        assert_eq!(
            GcContents::region_count_for_viewport(length, length - 1, None),
            2
        );
        assert_eq!(
            GcContents::region_count_for_viewport(length, length - 1, Some((length - 1, length))),
            1
        );
        assert_eq!(
            GcContents::region_count_for_viewport(length, 1, Some((length - 1, length))),
            1
        );
    }

    #[test]
    fn region_count_matches_materialized_bins_on_density_ladder_lengths() {
        for length in [20_000, 250_000, 2_000_000] {
            let sequence = b"aCGtN".repeat(length / 5);
            for bin_size in [100, 137, length + 1] {
                let gc = GcContents::new_from_sequence_with_bin_size(&sequence, bin_size);
                for viewport in [
                    None,
                    Some((0, 5_000)),
                    Some((1, 5_001)),
                    Some((99, 100)),
                    Some((100, 101)),
                    Some((length - 1, usize::MAX)),
                ] {
                    let expected = gc
                        .regions()
                        .iter()
                        .filter(|region| {
                            viewport.is_none_or(|(start, end)| {
                                region.from() < end && region.to() > start
                            })
                        })
                        .count();
                    assert_eq!(
                        GcContents::region_count_for_viewport(length, bin_size, viewport),
                        expected,
                        "length={length}, bin={bin_size}, viewport={viewport:?}"
                    );
                }
            }
        }
    }
}
