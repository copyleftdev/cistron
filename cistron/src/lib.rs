#![forbid(unsafe_code)]
//! `cistron` — the located allele and the rules that decide when two spellings
//! are the same variant.
//!
//! Coordinates are **interbase**: a position is the count of reference bases
//! before the locus, so insertions and deletions are symmetric and the
//! 0-based/1-based ambiguity is converted away at the edges rather than carried
//! through the core.
//!
//! Three normalizations live here, one per downstream convention:
//! [`Variant::normalize`] is LEFT-aligned + blunt/parsimonious (bcftools norm,
//! validated against it byte-for-byte); [`Variant::normalize_right`] is the
//! 3'-most dual (HGVS); [`Variant::fully_justified`] is GA4GH-VRS expansion
//! (validated against vrs-python). The canonical `normalize` is transcribed from
//! `specs/Normalize.tla` and held to it by the oracle harness in `../harness`.

/// A DNA base. `repr(u8)` so it round-trips through a compact integer encoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    /// The base for encoding `i` (`0=A,1=C,2=G,3=T`), or `None` if out of range.
    pub const fn from_index(i: u8) -> Option<Base> {
        match i {
            0 => Some(Base::A),
            1 => Some(Base::C),
            2 => Some(Base::G),
            3 => Some(Base::T),
            _ => None,
        }
    }

    pub const fn index(self) -> u8 {
        self as u8
    }
}

/// An interbase coordinate: the number of reference bases before a locus.
///
/// The field is private so a position can only be built through a
/// convention-named constructor — you cannot accidentally feed a 1-based VCF
/// `POS` where an interbase start is wanted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interbase(usize);

impl Interbase {
    /// An interbase coordinate stated directly (bases before the locus).
    pub const fn new(bases_before: usize) -> Self {
        Interbase(bases_before)
    }

    /// A BED / UCSC 0-based half-open **start** is already an interbase point.
    pub const fn from_zero_based_start(start: usize) -> Self {
        Interbase(start)
    }

    /// A VCF 1-based inclusive `POS` naming the first affected base sits one
    /// interbase point earlier.
    pub const fn from_one_based(pos: usize) -> Self {
        Interbase(pos - 1)
    }

    pub const fn get(self) -> usize {
        self.0
    }
}

/// A located allele edit against a reference: replace the `del` reference bases
/// starting at `pos` with `ins`. A pure insertion has empty `del`; a pure
/// deletion has empty `ins`; an empty–empty edit denotes no change.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Variant {
    pub pos: Interbase,
    pub del: Vec<Base>,
    pub ins: Vec<Base>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The variant reaches past the end of the reference.
    OutOfBounds { end: usize, reference_len: usize },
    /// `del` is not the reference sequence at the locus it claims to delete.
    ReferenceMismatch { pos: usize },
}

/// The fully-justified expansion of a variant — the raw material the VRS state
/// decision (literal vs run-length expression) reads. Interbase `[start, end)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VrsExpansion {
    pub start: usize,
    pub end: usize,
    /// The expanded alternate sequence.
    pub alt: Vec<Base>,
    /// Length of the trimmed unit (the deleted or inserted seed).
    pub seed_length: usize,
    /// Both alleles survived trimming — a substitution/complex edit, never a
    /// repeat expansion.
    pub both_sides: bool,
}

impl Variant {
    /// The alternate haplotype this variant produces against `reference`.
    /// Assumes the variant is in bounds (true once [`Variant::check`] passes).
    pub fn denotation(&self, reference: &[Base]) -> Vec<Base> {
        let pos = self.pos.get();
        let end = pos + self.del.len();
        let mut out = Vec::with_capacity(reference.len() + self.ins.len());
        out.extend_from_slice(&reference[..pos]);
        out.extend_from_slice(&self.ins);
        out.extend_from_slice(&reference[end..]);
        out
    }

    /// Well-formed iff in bounds and `del` is exactly the reference at the
    /// locus — reference-agreement, the check that catches the largest class of
    /// real-world variant bugs.
    pub fn check(&self, reference: &[Base]) -> Result<(), Error> {
        let pos = self.pos.get();
        let end = pos + self.del.len();
        if end > reference.len() {
            return Err(Error::OutOfBounds {
                end,
                reference_len: reference.len(),
            });
        }
        if self.del != reference[pos..end] {
            return Err(Error::ReferenceMismatch { pos });
        }
        Ok(())
    }

    /// True iff `del` and `ins` share no leading and no trailing base. A shared
    /// base is only possible when both alleles are non-empty, so any empty
    /// allele is trivially parsimonious.
    pub fn is_parsimonious(&self) -> bool {
        match (
            self.del.first(),
            self.del.last(),
            self.ins.first(),
            self.ins.last(),
        ) {
            (Some(df), Some(dl), Some(inf), Some(inl)) => df != inf && dl != inl,
            _ => true,
        }
    }

    /// The canonical form: LEFT-aligned + blunt/parsimonious. Idempotent, and
    /// denotation-preserving. Errors only if the variant fails [`Variant::check`].
    pub fn normalize(&self, reference: &[Base]) -> Result<Variant, Error> {
        self.check(reference)?;

        let mut pos = self.pos.get();
        let mut del = self.del.clone();
        let mut ins = self.ins.clone();

        // Trim the shared suffix, then the shared prefix (blunt form).
        while matches!((del.last(), ins.last()), (Some(a), Some(b)) if a == b) {
            del.pop();
            ins.pop();
        }
        while !del.is_empty() && !ins.is_empty() && del[0] == ins[0] {
            del.remove(0);
            ins.remove(0);
            pos += 1;
        }

        // Left-align: roll one base at a time while an equivalent lower-position
        // representation exists. The condition is local to the trimmed form.
        loop {
            if pos == 0 {
                break;
            }
            let dlen = del.len();
            let can_shift = if del.is_empty() && ins.is_empty() {
                true // null edit: only canonical at pos 0
            } else if ins.is_empty() {
                reference[pos - 1] == reference[pos - 1 + dlen] // pure deletion
            } else {
                *ins.last().unwrap() == reference[pos - 1 + dlen] // insertion / mixed
            };
            if !can_shift {
                break;
            }
            pos -= 1;
            if del.is_empty() && ins.is_empty() {
                continue;
            }
            let entering = reference[pos];
            if !del.is_empty() {
                del = reference[pos..pos + dlen].to_vec();
            }
            if !ins.is_empty() {
                ins.insert(0, entering);
                ins.pop();
            }
        }

        Ok(Variant {
            pos: Interbase(pos),
            del,
            ins,
        })
    }

    /// The **fully-justified** form (GA4GH VRS normalization): trim, then expand
    /// over the entire ambiguous region by rolling both left and right. Unlike
    /// [`Variant::normalize`] (leftmost, minimal), this widens the interval to
    /// cover a whole repeat, with the alternate sequence expanded to match — the
    /// representation VRS digests. Denotation-preserving.
    pub fn fully_justified(&self, reference: &[Base]) -> Result<Variant, Error> {
        let e = self.vrs_expand(reference)?;
        Ok(Variant {
            pos: Interbase(e.start),
            del: reference[e.start..e.end].to_vec(),
            ins: e.alt,
        })
    }

    /// The building blocks the VRS state decision (literal vs run-length) is
    /// derived from: the fully-justified interval and alt, plus the trimmed
    /// "seed" unit length and whether both alleles survived trimming (a
    /// substitution/complex edit, which VRS never expands to a repeat).
    pub fn vrs_expand(&self, reference: &[Base]) -> Result<VrsExpansion, Error> {
        self.check(reference)?;
        let mut start = self.pos.get();
        let mut end = start + self.del.len();
        let mut refa = self.del.clone();
        let mut alta = self.ins.clone();

        // Trim the common prefix, then the common suffix.
        let mut i = 0;
        while i < refa.len() && i < alta.len() && refa[i] == alta[i] {
            i += 1;
        }
        start += i;
        refa.drain(..i);
        alta.drain(..i);
        let mut j = 0;
        while j < refa.len()
            && j < alta.len()
            && refa[refa.len() - 1 - j] == alta[alta.len() - 1 - j]
        {
            j += 1;
        }
        end -= j;
        refa.truncate(refa.len() - j);
        alta.truncate(alta.len() - j);

        let seed_length = if refa.is_empty() {
            alta.len()
        } else {
            refa.len()
        };
        let both_sides = !refa.is_empty() && !alta.is_empty();

        // Expand: roll left and right as far as the alleles stay periodic with
        // the reference, then grow the interval and alt to cover the roll.
        let ldist = roll_left(reference, &refa, &alta, start);
        let rdist = roll_right(reference, &refa, &alta, end);
        let new_start = start - ldist;
        let new_end = end + rdist;
        let mut alt = Vec::with_capacity(ldist + alta.len() + rdist);
        alt.extend_from_slice(&reference[new_start..start]);
        alt.extend_from_slice(&alta);
        alt.extend_from_slice(&reference[end..new_end]);

        Ok(VrsExpansion {
            start: new_start,
            end: new_end,
            alt,
            seed_length,
            both_sides,
        })
    }

    /// The 3'-most blunt/parsimonious form — the mirror of [`Variant::normalize`],
    /// as HGVS requires. Same trimmed `(del, ins)` as the left form; only the
    /// position differs. Not the canonical id form; use [`Variant::normalize`]
    /// for that.
    pub fn normalize_right(&self, reference: &[Base]) -> Result<Variant, Error> {
        self.check(reference)?;
        let n = reference.len();

        let mut pos = self.pos.get();
        let mut del = self.del.clone();
        let mut ins = self.ins.clone();

        while matches!((del.last(), ins.last()), (Some(a), Some(b)) if a == b) {
            del.pop();
            ins.pop();
        }
        while !del.is_empty() && !ins.is_empty() && del[0] == ins[0] {
            del.remove(0);
            ins.remove(0);
            pos += 1;
        }

        loop {
            let dlen = del.len();
            if pos + dlen >= n {
                break;
            }
            let can_shift = if del.is_empty() && ins.is_empty() {
                true // null edit: only 3'-canonical at the sequence end
            } else if ins.is_empty() {
                reference[pos] == reference[pos + dlen] // pure deletion
            } else {
                *ins.first().unwrap() == reference[pos] // insertion / mixed
            };
            if !can_shift {
                break;
            }
            if del.is_empty() && ins.is_empty() {
                pos += 1;
                continue;
            }
            let entering = reference[pos + dlen];
            pos += 1;
            if !del.is_empty() {
                del = reference[pos..pos + dlen].to_vec();
            }
            if !ins.is_empty() {
                ins.remove(0);
                ins.push(entering);
            }
        }

        Ok(Variant {
            pos: Interbase(pos),
            del,
            ins,
        })
    }
}

/// Rolling distance to the left (bioutils EXPAND): how far both alleles stay
/// periodic with the reference reading leftward from `start`.
fn roll_left(seq: &[Base], ref_a: &[Base], alt_a: &[Base], start: usize) -> usize {
    if start == 0 {
        return 0;
    }
    let ref_pos = start - 1; // 0-based base just left of the interval
    let mut d = 0;
    while d <= ref_pos {
        let expected = seq[ref_pos - d];
        if !cycle_left_ok(ref_a, d, expected) || !cycle_left_ok(alt_a, d, expected) {
            break;
        }
        d += 1;
    }
    d
}

fn cycle_left_ok(a: &[Base], d: usize, expected: Base) -> bool {
    if a.is_empty() {
        return true;
    }
    let len = a.len() as isize;
    let idx = (-(d as isize + 1)).rem_euclid(len) as usize;
    a[idx] == expected
}

/// Rolling distance to the right, the mirror of [`roll_left`].
fn roll_right(seq: &[Base], ref_a: &[Base], alt_a: &[Base], end: usize) -> usize {
    let n = seq.len();
    if end >= n {
        return 0;
    }
    let max_d = (n - 1) - end;
    let mut d = 0;
    while d <= max_d {
        let expected = seq[end + d];
        if !cycle_right_ok(ref_a, d, expected) || !cycle_right_ok(alt_a, d, expected) {
            break;
        }
        d += 1;
    }
    d
}

fn cycle_right_ok(a: &[Base], d: usize, expected: Base) -> bool {
    if a.is_empty() {
        return true;
    }
    a[d % a.len()] == expected
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base(i: u8) -> Base {
        Base::from_index(i).unwrap()
    }

    /// Every base string of length `0..=max` over `{A, C}`.
    fn seqs(max: usize) -> Vec<Vec<Base>> {
        let mut out = vec![vec![]];
        let mut frontier = vec![vec![]];
        for _ in 0..max {
            let mut next = Vec::new();
            for s in &frontier {
                for b in [base(0), base(1), base(2)] {
                    let mut e = s.clone();
                    e.push(b);
                    next.push(e);
                }
            }
            out.extend(next.iter().cloned());
            frontier = next;
        }
        out
    }

    fn domain(max_ref: usize, max_allele: usize) -> Vec<(Vec<Base>, Variant)> {
        let mut out = Vec::new();
        for reference in seqs(max_ref).into_iter().filter(|r| !r.is_empty()) {
            for pos in 0..=reference.len() {
                let dmax = (reference.len() - pos).min(max_allele);
                for dlen in 0..=dmax {
                    let del = reference[pos..pos + dlen].to_vec();
                    for ins in seqs(max_allele) {
                        out.push((
                            reference.clone(),
                            Variant {
                                pos: Interbase::new(pos),
                                del: del.clone(),
                                ins,
                            },
                        ));
                    }
                }
            }
        }
        out
    }

    #[test]
    fn normalize_preserves_denotation_and_is_idempotent() {
        let mut checked = 0;
        for (reference, raw) in domain(4, 2) {
            let once = raw
                .normalize(&reference)
                .expect("well-formed input normalizes");
            assert_eq!(
                once.denotation(&reference),
                raw.denotation(&reference),
                "denotation changed for {raw:?}"
            );
            assert!(
                once.is_parsimonious(),
                "result not parsimonious for {raw:?}"
            );
            once.check(&reference).expect("result is well-formed");
            let twice = once.normalize(&reference).expect("re-normalizes");
            assert_eq!(twice, once, "not idempotent for {raw:?}");
            checked += 1;
        }
        assert!(checked > 100, "domain collapsed to {checked}");
    }

    #[test]
    fn base_encoding_round_trips() {
        for i in 0..4u8 {
            assert_eq!(Base::from_index(i).unwrap().index(), i);
        }
        assert_eq!(Base::from_index(0), Some(Base::A));
        assert_eq!(Base::from_index(3), Some(Base::T));
        assert_eq!(Base::from_index(4), None);
    }

    #[test]
    fn denotation_is_the_spliced_haplotype() {
        let reference = bases("ACGT");
        let snv = Variant {
            pos: Interbase::new(1),
            del: bases("C"),
            ins: bases("T"),
        };
        assert_eq!(snv.denotation(&reference), bases("ATGT"));
        let ins = Variant {
            pos: Interbase::new(2),
            del: vec![],
            ins: bases("A"),
        };
        assert_eq!(ins.denotation(&reference), bases("ACAGT"));
        let del = Variant {
            pos: Interbase::new(1),
            del: bases("C"),
            ins: vec![],
        };
        assert_eq!(del.denotation(&reference), bases("AGT"));
    }

    #[test]
    fn check_rejects_mismatch_and_out_of_bounds() {
        let reference = bases("ACGT");
        let good = Variant {
            pos: Interbase::new(1),
            del: bases("C"),
            ins: bases("T"),
        };
        assert_eq!(good.check(&reference), Ok(()));
        let wrong_ref = Variant {
            pos: Interbase::new(1),
            del: bases("A"),
            ins: bases("T"),
        };
        assert_eq!(
            wrong_ref.check(&reference),
            Err(Error::ReferenceMismatch { pos: 1 })
        );
        let past_end = Variant {
            pos: Interbase::new(3),
            del: bases("TT"),
            ins: vec![],
        };
        assert!(matches!(
            past_end.check(&reference),
            Err(Error::OutOfBounds { .. })
        ));
    }

    #[test]
    fn parsimony_detects_shared_prefix_and_suffix() {
        let shares_prefix = Variant {
            pos: Interbase::new(0),
            del: bases("CG"),
            ins: bases("CT"),
        };
        let shares_suffix = Variant {
            pos: Interbase::new(0),
            del: bases("GC"),
            ins: bases("TC"),
        };
        let distinct = Variant {
            pos: Interbase::new(0),
            del: bases("C"),
            ins: bases("T"),
        };
        let pure_del = Variant {
            pos: Interbase::new(0),
            del: bases("C"),
            ins: vec![],
        };
        let null = Variant {
            pos: Interbase::new(0),
            del: vec![],
            ins: vec![],
        };
        assert!(!shares_prefix.is_parsimonious());
        assert!(!shares_suffix.is_parsimonious());
        assert!(distinct.is_parsimonious());
        assert!(pure_del.is_parsimonious());
        assert!(null.is_parsimonious());
    }

    fn bases(s: &str) -> Vec<Base> {
        s.chars()
            .map(|c| match c {
                'A' => Base::A,
                'C' => Base::C,
                'G' => Base::G,
                'T' => Base::T,
                _ => panic!("bad base {c}"),
            })
            .collect()
    }

    #[test]
    fn interbase_conventions_agree() {
        // VCF POS=5 (1-based, first affected base) and BED start=4 name the
        // same interbase point.
        assert_eq!(
            Interbase::from_one_based(5),
            Interbase::from_zero_based_start(4)
        );
    }

    #[test]
    fn right_normalization_is_a_valid_3prime_dual() {
        let mut differ = 0;
        for (reference, raw) in domain(4, 2) {
            let left = raw.normalize(&reference).unwrap();
            let right = raw.normalize_right(&reference).unwrap();
            // The shared invariant is the denotation. Content may rotate in a
            // tandem repeat (AC<->CA), which is why the two conventions differ.
            assert_eq!(
                left.denotation(&reference),
                right.denotation(&reference),
                "denotation drifted for {raw:?}"
            );
            assert!(
                right.is_parsimonious(),
                "right form not parsimonious for {raw:?}"
            );
            assert!(
                left.pos.get() <= right.pos.get(),
                "left not <= right for {raw:?}"
            );
            assert_eq!(
                right.normalize_right(&reference).unwrap(),
                right,
                "right not idempotent for {raw:?}"
            );
            if left.pos != right.pos {
                differ += 1;
            }
        }
        assert!(differ > 0, "no variant shifted — repeats not exercised");
    }

    #[test]
    fn fully_justified_preserves_denotation_and_is_idempotent() {
        let mut expanded = 0;
        for (reference, raw) in domain(4, 2) {
            if raw.denotation(&reference) == reference {
                continue; // null: not a variant (vrs-python rejects these too)
            }
            let fj = raw.fully_justified(&reference).unwrap();
            assert_eq!(
                fj.denotation(&reference),
                raw.denotation(&reference),
                "denotation drifted for {raw:?}"
            );
            fj.check(&reference).expect("result is well-formed");
            assert_eq!(
                fj.fully_justified(&reference).unwrap(),
                fj,
                "not idempotent for {raw:?}"
            );
            if fj.del.len() > raw.del.len() {
                expanded += 1;
            }
        }
        assert!(expanded > 0, "no variant expanded — repeats not exercised");
    }

    #[test]
    fn vrs_expand_reports_seed_and_both_sides() {
        let reference = bases("ACGT");
        // substitution: both alleles survive trimming.
        let e = Variant {
            pos: Interbase::new(1),
            del: bases("C"),
            ins: bases("T"),
        }
        .vrs_expand(&reference)
        .unwrap();
        assert!(e.both_sides);
        assert_eq!(e.seed_length, 1);
        // pure deletion: not both sides; seed is the deleted length.
        let e = Variant {
            pos: Interbase::new(1),
            del: bases("C"),
            ins: vec![],
        }
        .vrs_expand(&reference)
        .unwrap();
        assert!(!e.both_sides);
        assert_eq!(e.seed_length, 1);
        // pure insertion: not both sides; seed is the inserted length.
        let e = Variant {
            pos: Interbase::new(1),
            del: vec![],
            ins: bases("CG"),
        }
        .vrs_expand(&reference)
        .unwrap();
        assert!(!e.both_sides);
        assert_eq!(e.seed_length, 2);
        // anchored form that trims down to a pure deletion.
        let e = Variant {
            pos: Interbase::new(0),
            del: bases("AC"),
            ins: bases("A"),
        }
        .vrs_expand(&reference)
        .unwrap();
        assert!(!e.both_sides);
    }

    #[test]
    fn fully_justified_golden_vectors() {
        // (reference, raw start/end/alt) -> (norm start/end/alt), from vrs-python.
        let cases: &[(&str, usize, usize, &str, usize, usize, &str)] = &[
            ("GCAAAAT", 5, 5, "A", 2, 6, "AAAAA"), // insert at run's right end rolls left
            ("CAAAAG", 2, 3, "", 1, 5, "AAA"),     // delete a middle base of a run
            ("GACACACT", 1, 1, "AC", 1, 7, "ACACACAC"), // dinucleotide insertion
            ("GACACACT", 1, 3, "", 1, 7, "ACAC"),  // dinucleotide deletion
            ("ACGT", 1, 2, "T", 1, 2, "T"),        // SNV: unchanged
            ("ACGT", 1, 3, "TT", 1, 3, "TT"),      // MNV: unchanged
            ("GAAAT", 1, 4, "CC", 1, 4, "CC"),     // complex delins: unchanged
            ("CGAT", 1, 3, "G", 2, 3, ""),         // anchored: shared prefix trims
            ("CTAG", 1, 3, "A", 1, 2, ""),         // anchored: shared suffix trims
            ("GCAGCAGCAGT", 7, 10, "", 0, 10, "GCAGCAG"), // trinucleotide left+right roll
            ("GCAGCAGCAGT", 4, 4, "CAG", 0, 10, "GCAGCAGCAGCAG"), // trinucleotide insertion
        ];
        for &(seq, s, e, alt, es, ee, ealt) in cases {
            let reference = bases(seq);
            let raw = Variant {
                pos: Interbase::new(s),
                del: reference[s..e].to_vec(),
                ins: bases(alt),
            };
            let fj = raw.fully_justified(&reference).unwrap();
            let got_alt: String = fj
                .ins
                .iter()
                .map(|b| ['A', 'C', 'G', 'T'][b.index() as usize])
                .collect();
            assert_eq!(
                (fj.pos.get(), fj.pos.get() + fj.del.len(), got_alt.as_str()),
                (es, ee, ealt),
                "fully_justified {seq} {s}:{e} {alt:?}"
            );
        }
    }

    /// Long homopolymers and multi-base repeat units: the multi-step roll must
    /// carry an indel all the way to the 5' end (left) and 3' end (right),
    /// preserving denotation, on sequences far longer than the exhaustive domain.
    #[test]
    fn long_repeats_roll_to_the_boundaries() {
        let a = base(0);
        let c = base(1);

        // Homopolymer A^L: any A-indel left-aligns to 0 and right-aligns to the
        // end of the run.
        for l in 2..=16usize {
            let reference = vec![a; l];
            for k in 1..=3.min(l) {
                for start in 0..=(l - k) {
                    let v = Variant {
                        pos: Interbase::new(start),
                        del: vec![a; k],
                        ins: vec![],
                    };
                    let left = v.normalize(&reference).unwrap();
                    let right = v.normalize_right(&reference).unwrap();
                    assert_eq!(left.pos.get(), 0, "A^{l} del{k}@{start}: left != 0");
                    assert_eq!(
                        right.pos.get(),
                        l - k,
                        "A^{l} del{k}@{start}: right != {}",
                        l - k
                    );
                    assert_eq!(left.denotation(&reference), v.denotation(&reference));
                    assert_eq!(right.denotation(&reference), v.denotation(&reference));
                }
            }
            for k in 1..=3usize {
                for start in 0..=l {
                    let v = Variant {
                        pos: Interbase::new(start),
                        del: vec![],
                        ins: vec![a; k],
                    };
                    let left = v.normalize(&reference).unwrap();
                    let right = v.normalize_right(&reference).unwrap();
                    assert_eq!(left.pos.get(), 0, "A^{l} ins{k}@{start}: left != 0");
                    assert_eq!(right.pos.get(), l, "A^{l} ins{k}@{start}: right != {l}");
                    assert_eq!(left.denotation(&reference), v.denotation(&reference));
                }
            }
        }

        // (AC)^k tandem repeat: deleting one unit rolls across a multi-base motif.
        for k in 2..=8usize {
            let reference: Vec<Base> = (0..k).flat_map(|_| [a, c]).collect();
            for u in 0..k {
                let start = u * 2;
                let v = Variant {
                    pos: Interbase::new(start),
                    del: vec![a, c],
                    ins: vec![],
                };
                let left = v.normalize(&reference).unwrap();
                let right = v.normalize_right(&reference).unwrap();
                assert_eq!(left.pos.get(), 0, "(AC)^{k} del@{start}: left != 0");
                assert_eq!(
                    right.pos.get(),
                    (k - 1) * 2,
                    "(AC)^{k} del@{start}: right != end"
                );
                assert_eq!(left.denotation(&reference), v.denotation(&reference));
                assert_eq!(right.denotation(&reference), v.denotation(&reference));
            }
        }
    }
}
