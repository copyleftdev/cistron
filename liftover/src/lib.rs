#![forbid(unsafe_code)]
//! Coordinate liftover between two assemblies, with one governing discipline:
//! **fallible, never silent.** A coordinate that cannot map cleanly returns a
//! typed error — never a plausible-but-wrong number. That is the whole reason
//! this is a separate rule: liftover is the one coordinate operation that is
//! genuinely partial (positions fall in gaps) and can invert strand.
//!
//! Coordinates are **interbase** (0-based, half-open), matching the rest of
//! `cistron`. A [`Chain`] is a list of aligned [`Block`]s; the gaps between
//! blocks are exactly where a lift fails. This crate is the chain *algebra* —
//! parsing UCSC chain files is a boundary concern layered on top, not here.

/// Which strand of the target the source aligns to.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Strand {
    Plus,
    Minus,
}

/// An aligned run: `len` source bases at interbase `src` map, in order, to
/// `len` target bases at interbase `tgt` (in the chain's own frame).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Block {
    pub src: usize,
    pub tgt: usize,
    pub len: usize,
}

/// An interbase interval `[start, end)` on the source.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

/// A lifted interval on the target, carrying the strand it landed on.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Lifted {
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
}

/// Why a lift produced no coordinate. Never a silently-wrong answer.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LiftError {
    /// The interval's start lies in a gap — no aligned target.
    Unmapped,
    /// The interval is not covered by a single aligned block, so it has no one
    /// contiguous image. (A boundary crate may choose to split it; the core
    /// refuses to guess.)
    Split,
    /// A zero-length interval; there is nothing to lift.
    Empty,
}

/// Why a chain could not be built.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BuildError {
    /// A block with zero length.
    EmptyBlock,
    /// Two blocks overlap on the source or target.
    Overlap,
    /// A block reaches past the source or target size.
    OutOfBounds,
}

/// An alignment chain from a source assembly to a target assembly.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Chain {
    blocks: Vec<Block>, // sorted by src, disjoint on both axes
    strand: Strand,
    src_size: usize,
    tgt_size: usize,
}

impl Chain {
    /// Build a chain, validating that blocks are non-empty, in bounds, and
    /// non-overlapping. Blocks are sorted by source position.
    pub fn new(
        src_size: usize,
        tgt_size: usize,
        strand: Strand,
        mut blocks: Vec<Block>,
    ) -> Result<Chain, BuildError> {
        blocks.sort_by_key(|b| b.src);
        let mut last_src_end = 0;
        let mut tgt_ranges: Vec<(usize, usize)> = Vec::with_capacity(blocks.len());
        for b in &blocks {
            if b.len == 0 {
                return Err(BuildError::EmptyBlock);
            }
            if b.src + b.len > src_size || b.tgt + b.len > tgt_size {
                return Err(BuildError::OutOfBounds);
            }
            if b.src < last_src_end {
                return Err(BuildError::Overlap);
            }
            last_src_end = b.src + b.len;
            let (s, e) = (b.tgt, b.tgt + b.len);
            if tgt_ranges.iter().any(|&(rs, re)| s < re && rs < e) {
                return Err(BuildError::Overlap);
            }
            tgt_ranges.push((s, e));
        }
        Ok(Chain {
            blocks,
            strand,
            src_size,
            tgt_size,
        })
    }

    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Lift a source interval to the target, or say precisely why it cannot.
    ///
    /// A non-empty interval maps only if a single block covers it end to end.
    /// If its start sits inside a block but it runs past the block's end, that
    /// is a [`LiftError::Split`]; if its start is in a gap, [`LiftError::Unmapped`].
    pub fn lift(&self, iv: Interval) -> Result<Lifted, LiftError> {
        if iv.start >= iv.end {
            return Err(LiftError::Empty);
        }
        // Fully covered by one block: the only clean case.
        for b in &self.blocks {
            if b.src <= iv.start && iv.end <= b.src + b.len {
                let cs = b.tgt + (iv.start - b.src);
                let ce = b.tgt + (iv.end - b.src);
                return Ok(match self.strand {
                    Strand::Plus => Lifted {
                        start: cs,
                        end: ce,
                        strand: Strand::Plus,
                    },
                    // Minus: reflect the target interval within the contig.
                    Strand::Minus => Lifted {
                        start: self.tgt_size - ce,
                        end: self.tgt_size - cs,
                        strand: Strand::Minus,
                    },
                });
            }
        }
        // Start covered but running past the block edge: a split, not a guess.
        for b in &self.blocks {
            if b.src <= iv.start && iv.start < b.src + b.len {
                return Err(LiftError::Split);
            }
        }
        Err(LiftError::Unmapped)
    }

    /// The inverse chain (target back to source), for Plus chains. Composing
    /// `lift` with `invert().lift` is the identity on the mapped domain.
    pub fn invert_plus(&self) -> Option<Chain> {
        if self.strand != Strand::Plus {
            return None;
        }
        let blocks = self
            .blocks
            .iter()
            .map(|b| Block {
                src: b.tgt,
                tgt: b.src,
                len: b.len,
            })
            .collect();
        Chain::new(self.tgt_size, self.src_size, Strand::Plus, blocks).ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chain(strand: Strand, src: usize, tgt: usize, blocks: &[(usize, usize, usize)]) -> Chain {
        let bs = blocks
            .iter()
            .map(|&(s, t, l)| Block {
                src: s,
                tgt: t,
                len: l,
            })
            .collect();
        Chain::new(src, tgt, strand, bs).unwrap()
    }

    #[test]
    fn plus_lift_is_a_simple_offset() {
        // src [0,10) maps to tgt [100,110); a gap in src [10,15).
        let c = chain(Strand::Plus, 20, 200, &[(0, 100, 10), (15, 110, 5)]);
        assert_eq!(
            c.lift(Interval { start: 2, end: 5 }),
            Ok(Lifted {
                start: 102,
                end: 105,
                strand: Strand::Plus
            })
        );
        assert_eq!(
            c.lift(Interval { start: 16, end: 18 }),
            Ok(Lifted {
                start: 111,
                end: 113,
                strand: Strand::Plus
            })
        );
    }

    #[test]
    fn gaps_are_unmapped_and_spans_are_split_never_wrong() {
        let c = chain(Strand::Plus, 20, 200, &[(0, 100, 10), (15, 110, 5)]);
        // Start inside the gap [10,15): unmapped.
        assert_eq!(
            c.lift(Interval { start: 11, end: 13 }),
            Err(LiftError::Unmapped)
        );
        // Start in block 0 but running into the gap: split.
        assert_eq!(
            c.lift(Interval { start: 8, end: 12 }),
            Err(LiftError::Split)
        );
        // Start exactly at a block's (exclusive) end is in the gap, not a split.
        assert_eq!(
            c.lift(Interval { start: 10, end: 12 }),
            Err(LiftError::Unmapped)
        );
        // Zero-length: nothing to lift.
        assert_eq!(c.lift(Interval { start: 5, end: 5 }), Err(LiftError::Empty));
    }

    #[test]
    fn minus_reflects_the_target_and_preserves_length() {
        // Same block on the minus strand of a target of size 200.
        let c = chain(Strand::Minus, 20, 200, &[(0, 100, 10)]);
        let lifted = c.lift(Interval { start: 2, end: 5 }).unwrap();
        // chain-frame target [102,105) reflects to [200-105, 200-102) = [95,98).
        assert_eq!(
            lifted,
            Lifted {
                start: 95,
                end: 98,
                strand: Strand::Minus
            }
        );
        // Length is invariant under lifting.
        assert_eq!(lifted.end - lifted.start, 3);
    }

    #[test]
    fn plus_round_trips_through_the_inverse() {
        let c = chain(
            Strand::Plus,
            50,
            300,
            &[(0, 100, 12), (20, 130, 8), (40, 200, 5)],
        );
        let inv = c.invert_plus().unwrap();
        for &(bs, _, len) in &[(0usize, 0usize, 12usize), (20, 0, 8), (40, 0, 5)] {
            for start in bs..bs + len {
                for end in start + 1..=bs + len {
                    let iv = Interval { start, end };
                    let up = c.lift(iv).unwrap();
                    let back = inv
                        .lift(Interval {
                            start: up.start,
                            end: up.end,
                        })
                        .unwrap();
                    assert_eq!(
                        (back.start, back.end),
                        (iv.start, iv.end),
                        "round-trip {iv:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn overlapping_or_out_of_bounds_blocks_are_rejected() {
        assert_eq!(
            Chain::new(
                10,
                10,
                Strand::Plus,
                vec![
                    Block {
                        src: 0,
                        tgt: 0,
                        len: 6
                    },
                    Block {
                        src: 4,
                        tgt: 8,
                        len: 2
                    }
                ]
            ),
            Err(BuildError::Overlap)
        );
        assert_eq!(
            Chain::new(
                10,
                10,
                Strand::Plus,
                vec![Block {
                    src: 0,
                    tgt: 8,
                    len: 5
                }]
            ),
            Err(BuildError::OutOfBounds)
        );
        assert_eq!(
            Chain::new(
                10,
                10,
                Strand::Plus,
                vec![Block {
                    src: 0,
                    tgt: 0,
                    len: 0
                }]
            ),
            Err(BuildError::EmptyBlock)
        );
    }
}
