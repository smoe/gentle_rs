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
}
