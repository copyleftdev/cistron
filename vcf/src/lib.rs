#![forbid(unsafe_code)]
//! The VCF boundary.
//!
//! VCF is 1-based and writes indels with a shared **anchor base** (`REF=GAT
//! ALT=G` deletes `AT`). cistron is interbase and blunt. Converting in:
//! `pos = POS - 1`, `del = REF`, `ins = ALT` — and [`cistron::Variant::normalize`]
//! trims the anchor for free. Converting out is where the work is: a blunt
//! insertion or deletion cannot be written in VCF without re-introducing an
//! anchor base drawn from the reference, so [`to_vcf`] needs the reference and
//! can fail when there is no base to anchor on.
//!
//! Symbolic ALTs (`<DEL>`), breakends (`[`/`]`), the spanning-deletion `*`, and
//! any non-ACGT base are outside the literal-allele core and are reported as
//! errors, never silently mangled.

use cistron::{Base, Interbase, Variant};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VcfError {
    /// Fewer than the five mandatory columns (CHROM POS ID REF ALT).
    TooFewColumns,
    /// POS was not a positive integer.
    BadPos(String),
    /// REF or ALT was empty; VCF requires at least one base.
    EmptyAllele,
    /// A base outside {A,C,G,T} (e.g. `N`).
    NonAcgt(char),
    /// A symbolic, breakend, spanning, or missing ALT — not a literal allele.
    UnsupportedAlt(String),
    /// The variant denotes no change; VCF has no representation for it.
    NullVariant,
    /// A blunt indel at the very start of the reference with no flanking base
    /// to anchor on.
    CannotAnchor,
}

/// A cistron variant with the CHROM/ID context a VCF row carries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocatedVariant {
    pub chrom: String,
    pub id: String,
    pub variant: Variant,
}

/// A single VCF allele: a 1-based position and the anchored REF/ALT bases.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfAllele {
    pub pos: u64,
    pub reference: Vec<Base>,
    pub alt: Vec<Base>,
}

impl VcfAllele {
    pub fn ref_str(&self) -> String {
        render(&self.reference)
    }
    pub fn alt_str(&self) -> String {
        render(&self.alt)
    }
}

/// Parse one `(POS, REF, ALT)` literal allele into an interbase variant. The
/// result is *not* normalized — the anchor base is still present until you call
/// [`cistron::Variant::normalize`].
pub fn to_variant(pos: u64, reference_field: &str, alt_field: &str) -> Result<Variant, VcfError> {
    if pos == 0 {
        return Err(VcfError::BadPos("0".into()));
    }
    let del = parse_bases(reference_field)?;
    let ins = parse_alt(alt_field)?;
    if del.is_empty() {
        return Err(VcfError::EmptyAllele);
    }
    Ok(Variant {
        pos: Interbase::from_one_based(pos as usize),
        del,
        ins,
    })
}

/// Parse one VCF data line, splitting a multiallelic ALT into one variant per
/// alternate allele. Symbolic/breakend/non-ACGT alternates surface as errors.
pub fn parse_line(line: &str) -> Result<Vec<LocatedVariant>, VcfError> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 5 {
        return Err(VcfError::TooFewColumns);
    }
    let chrom = cols[0].to_string();
    let pos: u64 = cols[1]
        .parse()
        .map_err(|_| VcfError::BadPos(cols[1].into()))?;
    let id = cols[2].to_string();
    let reference_field = cols[3];

    let mut out = Vec::new();
    for alt in cols[4].split(',') {
        let variant = to_variant(pos, reference_field, alt)?;
        out.push(LocatedVariant {
            chrom: chrom.clone(),
            id: id.clone(),
            variant,
        });
    }
    Ok(out)
}

/// Emit a variant as a VCF allele against `reference`, re-introducing an anchor
/// base where a blunt indel needs one. Prefers a left anchor; falls back to a
/// right anchor only at the start of the reference.
pub fn to_vcf(reference: &[Base], v: &Variant) -> Result<VcfAllele, VcfError> {
    let pos = v.pos.get();

    if !v.del.is_empty() && !v.ins.is_empty() {
        // A substitution or complex indel: parsimonious, so no anchor needed.
        return Ok(VcfAllele {
            pos: (pos + 1) as u64,
            reference: v.del.clone(),
            alt: v.ins.clone(),
        });
    }
    if v.del.is_empty() && v.ins.is_empty() {
        return Err(VcfError::NullVariant);
    }

    if pos >= 1 {
        let anchor = reference[pos - 1];
        Ok(VcfAllele {
            pos: pos as u64,
            reference: prepend(anchor, &v.del),
            alt: prepend(anchor, &v.ins),
        })
    } else {
        // At the start: the anchor must come from after the locus.
        let after = v.del.len();
        if after >= reference.len() {
            return Err(VcfError::CannotAnchor);
        }
        let anchor = reference[after];
        Ok(VcfAllele {
            pos: 1,
            reference: append(&v.del, anchor),
            alt: append(&v.ins, anchor),
        })
    }
}

fn prepend(b: Base, rest: &[Base]) -> Vec<Base> {
    let mut v = Vec::with_capacity(rest.len() + 1);
    v.push(b);
    v.extend_from_slice(rest);
    v
}

fn append(rest: &[Base], b: Base) -> Vec<Base> {
    let mut v = rest.to_vec();
    v.push(b);
    v
}

fn parse_bases(s: &str) -> Result<Vec<Base>, VcfError> {
    s.chars()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => Ok(Base::A),
            'C' => Ok(Base::C),
            'G' => Ok(Base::G),
            'T' => Ok(Base::T),
            other => Err(VcfError::NonAcgt(other)),
        })
        .collect()
}

fn parse_alt(alt: &str) -> Result<Vec<Base>, VcfError> {
    if alt.is_empty() {
        return Err(VcfError::EmptyAllele);
    }
    if alt == "." || alt == "*" || alt.starts_with('<') || alt.contains('[') || alt.contains(']') {
        return Err(VcfError::UnsupportedAlt(alt.into()));
    }
    parse_bases(alt)
}

fn render(bases: &[Base]) -> String {
    bases
        .iter()
        .map(|b| match b {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base(i: u8) -> Base {
        Base::from_index(i).unwrap()
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

    /// The round-trip law: emitting a normalized variant to VCF and parsing it
    /// back lands on the same normal form. VCF is lossless up to normalization.
    #[test]
    fn emit_then_parse_preserves_the_normal_form() {
        let mut round_tripped = 0usize;
        let mut null = 0usize;
        let mut unanchorable = 0usize;

        for (reference, raw) in domain(3, 2) {
            let norm = raw.normalize(&reference).unwrap();
            let allele = match to_vcf(&reference, &norm) {
                Ok(a) => a,
                Err(VcfError::NullVariant) => {
                    null += 1;
                    continue;
                }
                Err(VcfError::CannotAnchor) => {
                    unanchorable += 1;
                    continue;
                }
                Err(e) => panic!("unexpected emit error {e:?} for {norm:?}"),
            };

            // Parse back through the text path to exercise the real parser.
            let reparsed = to_variant(allele.pos, &allele.ref_str(), &allele.alt_str())
                .expect("emitted allele parses");
            assert_eq!(
                reparsed.denotation(&reference),
                norm.denotation(&reference),
                "denotation drifted for {norm:?} via {allele:?}"
            );
            assert_eq!(
                reparsed.normalize(&reference).unwrap(),
                norm,
                "normal form drifted for {norm:?} via {allele:?}"
            );
            round_tripped += 1;
        }

        eprintln!(
            "round-tripped {round_tripped} variants through VCF \
             ({null} null, {unanchorable} unanchorable at reference start)"
        );
        assert!(round_tripped > 100, "domain collapsed to {round_tripped}");
    }

    #[test]
    fn parse_splits_multiallelic() {
        // chrom pos id REF=A ALT=T,C
        let line = "chr1\t5\trs1\tA\tT,C";
        let variants = parse_line(line).unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].variant.ins, vec![Base::T]);
        assert_eq!(variants[1].variant.ins, vec![Base::C]);
        assert_eq!(variants[0].variant.pos, Interbase::from_one_based(5));
        assert_eq!(variants[0].chrom, "chr1");
    }

    #[test]
    fn deletion_anchor_normalizes_away() {
        // POS=3 REF=GAT ALT=G on ref ...GAT... deletes AT after G.
        let reference = vec![Base::A, Base::A, Base::G, Base::A, Base::T]; // AAGAT
        let v = to_variant(3, "GAT", "G").unwrap();
        assert_eq!(v.pos, Interbase::from_one_based(3)); // interbase 2
        let norm = v.normalize(&reference).unwrap();
        // anchor G trimmed -> deletes AT at interbase 3
        assert_eq!(norm.ins, Vec::<Base>::new());
        assert_eq!(norm.del, vec![Base::A, Base::T]);
        assert_eq!(norm.pos, Interbase::new(3));
    }

    #[test]
    fn symbolic_and_breakend_and_n_are_rejected() {
        assert!(matches!(
            to_variant(10, "A", "<DEL>"),
            Err(VcfError::UnsupportedAlt(_))
        ));
        assert!(matches!(
            to_variant(10, "A", "A[chr2:20["),
            Err(VcfError::UnsupportedAlt(_))
        ));
        assert!(matches!(
            to_variant(10, "A", "*"),
            Err(VcfError::UnsupportedAlt(_))
        ));
        assert!(matches!(
            to_variant(10, "N", "A"),
            Err(VcfError::NonAcgt('N'))
        ));
        assert!(matches!(
            to_variant(10, "A", "."),
            Err(VcfError::UnsupportedAlt(_))
        ));
    }

    #[test]
    fn too_few_columns_is_an_error() {
        assert!(matches!(
            parse_line("chr1\t5\trs1\tA"),
            Err(VcfError::TooFewColumns)
        ));
    }

    #[test]
    fn empty_alt_is_rejected() {
        assert!(matches!(to_variant(5, "A", ""), Err(VcfError::EmptyAllele)));
    }

    #[test]
    fn to_vcf_left_anchors_a_pure_deletion() {
        // ACGT, delete C at interbase 1 -> left anchor A -> POS=1 REF=AC ALT=A.
        let reference = vec![Base::A, Base::C, Base::G, Base::T];
        let del = Variant {
            pos: Interbase::new(1),
            del: vec![Base::C],
            ins: vec![],
        };
        let a = to_vcf(&reference, &del).unwrap();
        assert_eq!(
            (a.pos, a.ref_str(), a.alt_str()),
            (1, "AC".into(), "A".into())
        );
    }
}
