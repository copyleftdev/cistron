#![forbid(unsafe_code)]
//! The HGVS genomic (`g.`) boundary.
//!
//! HGVS is the mirror of the VCF/VRS side in two ways that this crate has to
//! honor. It normalizes to the **3'-most** position (right-shift, via
//! [`cistron::Variant::normalize_right`]) where VCF left-aligns; and it **omits
//! the deleted/duplicated bases** from the string, so parsing a `del`/`dup`
//! needs the reference to recover them. It also distinguishes a tandem
//! duplication (`dup`) from a plain insertion (`ins`) — a rule that only makes
//! sense once the variant is 3'-shifted.
//!
//! Scope is the genomic level only, covering substitutions, deletions,
//! insertions, duplications, delins, and inversions. Transcript/protein levels
//! (`c.`/`p.`/`n.`/`m.`) and conversions are out of the core and are rejected,
//! not guessed. Genomic rendering is validated against biocommons `hgvs`.

use cistron::{Base, Interbase, Variant};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HgvsError {
    /// Not a genomic (`g.`) expression.
    NotGenomic(String),
    /// The expression did not parse as a supported operation.
    BadSyntax(String),
    /// A recognized but unsupported operation (inversion, conversion, ...).
    Unsupported(String),
    /// A base outside {A,C,G,T}.
    NonAcgt(char),
    /// A coordinate outside the reference.
    OutOfBounds,
    /// The stated reference base disagrees with the reference.
    ReferenceMismatch,
    /// `g.=` — no change; not a variant.
    NoChange,
    /// A blunt insertion at interbase 0 has no 5' flank to name in `g.`.
    AtSequenceStart,
}

/// Render a 3'-normalized variant as an HGVS genomic expression.
pub fn to_hgvs(reference: &[Base], variant: &Variant) -> Result<String, HgvsError> {
    let r = variant
        .normalize_right(reference)
        .map_err(|_| HgvsError::ReferenceMismatch)?;
    let pos = r.pos.get();
    let d = r.del.len();
    let m = r.ins.len();

    if d == 0 && m == 0 {
        return Err(HgvsError::NoChange);
    }

    let start1 = pos + 1;
    let end1 = pos + d;

    if d == 1 && m == 1 {
        return Ok(format!("g.{start1}{}>{}", ch(r.del[0]), ch(r.ins[0])));
    }
    if m == 0 {
        return Ok(if d == 1 {
            format!("g.{start1}del")
        } else {
            format!("g.{start1}_{end1}del")
        });
    }
    if d == 0 {
        // Duplication: the m bases immediately 5' of the (3'-shifted) insertion
        // point are exactly what is being inserted.
        if pos >= m && reference[pos - m..pos] == r.ins[..] {
            let dup_end = pos;
            let dup_start = pos - m + 1;
            return Ok(if m == 1 {
                format!("g.{dup_end}dup")
            } else {
                format!("g.{dup_start}_{dup_end}dup")
            });
        }
        if pos == 0 {
            // No 5' base to flank an insertion; HGVS expresses it as a delins on
            // the first reference base (insert the bases before it).
            return match reference.first() {
                Some(&anchor) => Ok(format!("g.1delins{}{}", render(&r.ins), ch(anchor))),
                None => Err(HgvsError::AtSequenceStart),
            };
        }
        return Ok(format!("g.{pos}_{}ins{}", pos + 1, render(&r.ins)));
    }
    // Inversion: the alternate is the reverse-complement of the deleted bases
    // (a block replaced by its own reverse complement). HGVS prefers `inv`.
    if d >= 2 && is_reverse_complement(&r.del, &r.ins) {
        return Ok(format!("g.{start1}_{end1}inv"));
    }
    Ok(if d == 1 {
        format!("g.{start1}delins{}", render(&r.ins))
    } else {
        format!("g.{start1}_{end1}delins{}", render(&r.ins))
    })
}

fn complement(b: Base) -> Base {
    match b {
        Base::A => Base::T,
        Base::T => Base::A,
        Base::C => Base::G,
        Base::G => Base::C,
    }
}

fn is_reverse_complement(del: &[Base], ins: &[Base]) -> bool {
    del.len() == ins.len()
        && del
            .iter()
            .rev()
            .map(|&b| complement(b))
            .eq(ins.iter().copied())
}

/// Parse an HGVS genomic expression into an interbase variant. Needs the
/// reference because HGVS omits deleted and duplicated bases.
pub fn from_hgvs(reference: &[Base], expr: &str) -> Result<Variant, HgvsError> {
    let body = expr
        .strip_prefix("g.")
        .ok_or_else(|| HgvsError::NotGenomic(expr.into()))?;
    let n = reference.len();

    let (a, rest) = take_uint(body).ok_or_else(|| HgvsError::BadSyntax(body.into()))?;
    let (b, tail, ranged) = match rest.strip_prefix('_') {
        Some(r) => {
            let (b, t) = take_uint(r).ok_or_else(|| HgvsError::BadSyntax(body.into()))?;
            (b, t, true)
        }
        None => (a, rest, false),
    };

    if tail.contains('>') {
        if ranged {
            return Err(HgvsError::BadSyntax(body.into()));
        }
        let cs: Vec<char> = tail.chars().collect();
        if cs.len() != 3 || cs[1] != '>' {
            return Err(HgvsError::BadSyntax(body.into()));
        }
        let (refb, altb) = (base_of(cs[0])?, base_of(cs[2])?);
        if a == 0 || a > n {
            return Err(HgvsError::OutOfBounds);
        }
        if reference[a - 1] != refb {
            return Err(HgvsError::ReferenceMismatch);
        }
        return Ok(Variant {
            pos: Interbase::from_one_based(a),
            del: vec![refb],
            ins: vec![altb],
        });
    }
    if let Some(seq) = tail.strip_prefix("delins") {
        let del = ref_slice(reference, a, b, n)?;
        return Ok(Variant {
            pos: Interbase::from_one_based(a),
            del,
            ins: parse_seq(seq)?,
        });
    }
    if tail == "del" {
        let del = ref_slice(reference, a, b, n)?;
        return Ok(Variant {
            pos: Interbase::from_one_based(a),
            del,
            ins: vec![],
        });
    }
    if let Some(seq) = tail.strip_prefix("ins") {
        // `a.checked_add(1)` so a coordinate near usize::MAX can't overflow.
        if !ranged || a > n || a.checked_add(1) != Some(b) {
            return Err(HgvsError::BadSyntax(body.into()));
        }
        let ins = parse_seq(seq)?;
        if ins.is_empty() {
            return Err(HgvsError::BadSyntax(body.into()));
        }
        return Ok(Variant {
            pos: Interbase::new(a),
            del: vec![],
            ins,
        });
    }
    if tail == "dup" {
        let ins = ref_slice(reference, a, b, n)?;
        return Ok(Variant {
            pos: Interbase::new(b),
            del: vec![],
            ins,
        });
    }
    if tail == "inv" {
        // The segment is replaced by its reverse complement.
        let del = ref_slice(reference, a, b, n)?;
        let ins = del.iter().rev().map(|&x| complement(x)).collect();
        return Ok(Variant {
            pos: Interbase::from_one_based(a),
            del,
            ins,
        });
    }
    if tail == "=" {
        return Err(HgvsError::NoChange);
    }
    Err(HgvsError::Unsupported(body.into()))
}

/// The reference bases at 1-based inclusive `[a, b]`, bounds-checked.
fn ref_slice(reference: &[Base], a: usize, b: usize, n: usize) -> Result<Vec<Base>, HgvsError> {
    if a == 0 || a > b || b > n {
        return Err(HgvsError::OutOfBounds);
    }
    Ok(reference[a - 1..b].to_vec())
}

fn take_uint(s: &str) -> Option<(usize, &str)> {
    let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
    if end == 0 {
        return None;
    }
    Some((s[..end].parse().ok()?, &s[end..]))
}

fn base_of(c: char) -> Result<Base, HgvsError> {
    match c.to_ascii_uppercase() {
        'A' => Ok(Base::A),
        'C' => Ok(Base::C),
        'G' => Ok(Base::G),
        'T' => Ok(Base::T),
        other => Err(HgvsError::NonAcgt(other)),
    }
}

fn parse_seq(s: &str) -> Result<Vec<Base>, HgvsError> {
    s.chars().map(base_of).collect()
}

fn ch(b: Base) -> char {
    match b {
        Base::A => 'A',
        Base::C => 'C',
        Base::G => 'G',
        Base::T => 'T',
    }
}

fn render(bases: &[Base]) -> String {
    bases.iter().map(|&b| ch(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base(i: u8) -> Base {
        Base::from_index(i).unwrap()
    }

    fn bases(s: &str) -> Vec<Base> {
        s.chars().map(|c| base_of(c).unwrap()).collect()
    }

    fn seqs(max: usize) -> Vec<Vec<Base>> {
        let mut out = vec![vec![]];
        let mut frontier = vec![vec![]];
        for _ in 0..max {
            let mut next = Vec::new();
            for s in &frontier {
                for b in [base(0), base(1)] {
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

    /// Round-trip: render to HGVS, parse back, and the LEFT canonical is
    /// recovered. Proves the 3'-shift and the dup/ins rule are self-consistent
    /// with the rest of cistron.
    #[test]
    fn hgvs_round_trips_to_the_canonical_form() {
        let mut round_tripped = 0;
        let mut no_change = 0;
        let mut at_start = 0;

        for (reference, raw) in domain(3, 2) {
            let left = raw.normalize(&reference).unwrap();
            let expr = match to_hgvs(&reference, &raw) {
                Ok(e) => e,
                Err(HgvsError::NoChange) => {
                    no_change += 1;
                    continue;
                }
                Err(HgvsError::AtSequenceStart) => {
                    at_start += 1;
                    continue;
                }
                Err(e) => panic!("unexpected emit error {e:?} for {raw:?}"),
            };
            let parsed = from_hgvs(&reference, &expr)
                .unwrap_or_else(|e| panic!("failed to parse {expr:?}: {e:?}"));
            assert_eq!(
                parsed.denotation(&reference),
                left.denotation(&reference),
                "denotation drifted: {raw:?} -> {expr}"
            );
            assert_eq!(
                parsed.normalize(&reference).unwrap(),
                left,
                "did not recover canonical: {raw:?} -> {expr}"
            );
            round_tripped += 1;
        }

        eprintln!(
            "round-tripped {round_tripped} variants through HGVS \
             ({no_change} no-change, {at_start} start-insertions)"
        );
        assert!(round_tripped > 100, "domain collapsed to {round_tripped}");
    }

    #[test]
    fn deletion_is_three_prime_shifted() {
        // CAAAG: deleting an A is reported at the 3'-most A (position 4).
        let reference = bases("CAAAG");
        let del_first_a = Variant {
            pos: Interbase::from_one_based(2),
            del: bases("A"),
            ins: vec![],
        };
        assert_eq!(to_hgvs(&reference, &del_first_a).unwrap(), "g.4del");
    }

    #[test]
    fn tandem_insertion_is_a_dup() {
        // CAAG: inserting an A into the AA run is a duplication of the 3'-most A.
        let reference = bases("CAAG");
        let insert_a = Variant {
            pos: Interbase::new(1),
            del: vec![],
            ins: bases("A"),
        };
        assert_eq!(to_hgvs(&reference, &insert_a).unwrap(), "g.3dup");
    }

    #[test]
    fn substitution_and_delins_render() {
        let reference = bases("ACGT");
        let snv = Variant {
            pos: Interbase::from_one_based(2),
            del: bases("C"),
            ins: bases("T"),
        };
        assert_eq!(to_hgvs(&reference, &snv).unwrap(), "g.2C>T");
        let delins = Variant {
            pos: Interbase::from_one_based(2),
            del: bases("CG"),
            ins: bases("A"),
        };
        assert_eq!(to_hgvs(&reference, &delins).unwrap(), "g.2_3delinsA");
    }

    /// Regression for a fuzz-found overflow: a coordinate near usize::MAX must
    /// not panic when the parser checks `b == a + 1`.
    #[test]
    fn parser_survives_overflow_coordinates() {
        let reference = bases("ACGT");
        assert!(from_hgvs(&reference, &format!("g.{0}_{0}insA", usize::MAX)).is_err());
        assert!(from_hgvs(&reference, &format!("g.{}del", usize::MAX)).is_err());
        assert!(from_hgvs(&reference, &format!("g.{}A>T", usize::MAX)).is_err());
    }

    #[test]
    fn non_genomic_and_unsupported_are_rejected() {
        let reference = bases("ACGT");
        assert!(matches!(
            from_hgvs(&reference, "c.2C>T"),
            Err(HgvsError::NotGenomic(_))
        ));
        // conversions remain out of scope
        assert!(matches!(
            from_hgvs(&reference, "g.1_4con"),
            Err(HgvsError::Unsupported(_))
        ));
        assert!(matches!(
            from_hgvs(&reference, "g.2C>N"),
            Err(HgvsError::NonAcgt('N'))
        ));
        assert!(matches!(
            from_hgvs(&reference, "g.9del"),
            Err(HgvsError::OutOfBounds)
        ));
    }

    #[test]
    fn inversion_round_trips() {
        // ACGT[1,4) inverted: ref "ACG" -> revcomp "CGT".
        let reference = bases("ACGT");
        let inv = Variant {
            pos: Interbase::new(0),
            del: bases("ACG"),
            ins: bases("CGT"),
        };
        assert_eq!(to_hgvs(&reference, &inv).unwrap(), "g.1_3inv");
        let parsed = from_hgvs(&reference, "g.1_3inv").unwrap();
        assert_eq!(parsed.denotation(&reference), inv.denotation(&reference));
        assert_eq!(parsed, inv);
    }

    #[test]
    fn insertion_flanks_must_be_consecutive() {
        let reference = bases("ACGT");
        // non-consecutive flanks
        assert!(matches!(
            from_hgvs(&reference, "g.2_5insA"),
            Err(HgvsError::BadSyntax(_))
        ));
        // an insertion needs a range at all
        assert!(matches!(
            from_hgvs(&reference, "g.2insA"),
            Err(HgvsError::BadSyntax(_))
        ));
        // consecutive but past the end of the reference
        assert!(matches!(
            from_hgvs(&reference, "g.5_6insA"),
            Err(HgvsError::BadSyntax(_))
        ));
    }

    #[test]
    fn substitution_must_be_exactly_ref_gt_alt() {
        let reference = bases("ACGT");
        assert!(matches!(
            from_hgvs(&reference, "g.2A>TT"),
            Err(HgvsError::BadSyntax(_))
        ));
    }

    #[test]
    fn to_hgvs_renders_a_plain_insertion() {
        // AT, insert C between them: not a duplication, not at the start.
        let reference = bases("AT");
        let ins = Variant {
            pos: Interbase::new(1),
            del: vec![],
            ins: bases("C"),
        };
        assert_eq!(to_hgvs(&reference, &ins).unwrap(), "g.1_2insC");
    }

    #[test]
    fn range_bounds_are_checked() {
        let reference = bases("ACGT");
        assert!(matches!(
            from_hgvs(&reference, "g.3_2del"), // start > end
            Err(HgvsError::OutOfBounds)
        ));
        assert!(matches!(
            from_hgvs(&reference, "g.3_9del"), // end past reference
            Err(HgvsError::OutOfBounds)
        ));
        assert!(matches!(
            from_hgvs(&reference, "g.0del"), // start at 0
            Err(HgvsError::OutOfBounds)
        ));
    }
}
