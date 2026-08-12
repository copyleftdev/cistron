#![forbid(unsafe_code)]
//! Content-addressed variant identifiers.
//!
//! A variant's identity is the digest of its *normalized* form, so two
//! spellings of the same edit hash to the same id and equality becomes a hash
//! compare. This follows the GA4GH-VRS shape: a variant is a location on a
//! named sequence plus the alternate state, and the identifier digests exactly
//! that. Following VRS, the deleted bases are **not** hashed — after
//! normalization they are fully determined by `(sequence, start, end)`.
//!
//! Two id schemes live here. [`variant_id`] is a compact internal
//! content-address (`cistron:va.`) over `cistron`'s own canonical bytes — useful
//! as a dedup/join key. The [`vrs`] module and [`ga4gh_allele_id`] emit real
//! **GA4GH VRS** identifiers (`ga4gh:VA.`), byte-for-byte with the spec (proven
//! against its worked example), so `cistron` interoperates with every other VRS
//! tool.

pub mod vrs;
use cistron::{Error, Variant};
use sha2::{Digest, Sha512};

/// A stable, content-addressed identifier for a normalized variant.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VariantId(String);

impl VariantId {
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::fmt::Display for VariantId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// The bytes that get digested: a canonical, length-delimited encoding of the
/// normalized variant on its named sequence. `sequence_id` is the namespace
/// (a RefSeq accession, a VRS `SQ.` digest, ...) that pins the reference —
/// without it the same coordinates on two assemblies would collide.
pub fn canonical_bytes(
    sequence_id: &str,
    reference: &[cistron::Base],
    variant: &Variant,
) -> Result<Vec<u8>, Error> {
    let norm = variant.normalize(reference)?;
    let start = norm.pos.get() as u64;
    let end = start + norm.del.len() as u64;

    let mut buf = Vec::new();
    put_bytes(&mut buf, sequence_id.as_bytes());
    buf.extend_from_slice(&start.to_be_bytes());
    buf.extend_from_slice(&end.to_be_bytes());
    let alt: Vec<u8> = norm.ins.iter().map(|b| b.index()).collect();
    put_bytes(&mut buf, &alt);
    Ok(buf)
}

/// The content-addressed id of `variant` on the sequence named `sequence_id`.
pub fn variant_id(
    sequence_id: &str,
    reference: &[cistron::Base],
    variant: &Variant,
) -> Result<VariantId, Error> {
    let bytes = canonical_bytes(sequence_id, reference, variant)?;
    let digest = Sha512::digest(&bytes);
    Ok(VariantId(format!(
        "cistron:va.{}",
        base64url_nopad(&digest[..24])
    )))
}

/// The real GA4GH VRS Allele identifier (`ga4gh:VA.…`) for `variant` on the
/// sequence named by its refget accession (`refget_accession`, carrying its own
/// `SQ.` prefix). The variant is left-aligned by [`cistron::Variant::normalize`]
/// first.
///
/// This matches vrs-python for every variant whose left-aligned form equals its
/// VRS "fully justified" form — all substitutions and non-repeat indels. Indels
/// inside a tandem repeat need VRS fully-justified normalization to agree; for
/// those, pass VRS-normalized `start`/`end`/`alt` straight to [`vrs::allele_id`]
/// until that normalization lands.
pub fn ga4gh_allele_id(
    refget_accession: &str,
    reference: &[cistron::Base],
    variant: &Variant,
) -> Result<String, Error> {
    let norm = variant.normalize(reference)?;
    let start = norm.pos.get() as u64;
    let end = start + norm.del.len() as u64;
    let alt: String = norm
        .ins
        .iter()
        .map(|b| ['A', 'C', 'G', 'T'][b.index() as usize])
        .collect();
    Ok(vrs::allele_id(refget_accession, start, end, &alt))
}

/// Length-prefixed (u32 big-endian) so no two field boundaries can be confused.
fn put_bytes(buf: &mut Vec<u8>, bytes: &[u8]) {
    buf.extend_from_slice(&(bytes.len() as u32).to_be_bytes());
    buf.extend_from_slice(bytes);
}

pub(crate) fn base64url_nopad(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        // Pack up to three bytes big-endian into the low 24 bits.
        let n = u32::from_be_bytes([
            0,
            chunk[0],
            *chunk.get(1).unwrap_or(&0),
            *chunk.get(2).unwrap_or(&0),
        ]);
        out.push(ALPHABET[(n >> 18 & 63) as usize] as char);
        out.push(ALPHABET[(n >> 12 & 63) as usize] as char);
        if chunk.len() > 1 {
            out.push(ALPHABET[(n >> 6 & 63) as usize] as char);
        }
        if chunk.len() > 2 {
            out.push(ALPHABET[(n & 63) as usize] as char);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use cistron::{Base, Interbase, Variant};

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

    /// The identity law: on one sequence, two variants share an id **iff** they
    /// denote the same edit. Equality by hash is therefore sound (no false
    /// merges) and injective over the domain (no accidental collisions).
    #[test]
    fn id_equality_matches_denotation() {
        let seq = "cistron:test-contig";
        for (reference, variants) in by_reference(domain(3, 2)) {
            let entries: Vec<(VariantId, Vec<Base>)> = variants
                .iter()
                .map(|v| {
                    (
                        variant_id(seq, &reference, v).unwrap(),
                        v.denotation(&reference),
                    )
                })
                .collect();
            for (id_a, den_a) in &entries {
                for (id_b, den_b) in &entries {
                    assert_eq!(
                        id_a == id_b,
                        den_a == den_b,
                        "id/denotation disagree on ref={reference:?}: {id_a} vs {id_b}"
                    );
                }
            }
        }
    }

    /// Group `(reference, variant)` pairs by their reference.
    fn by_reference(pairs: Vec<(Vec<Base>, Variant)>) -> Vec<(Vec<Base>, Vec<Variant>)> {
        let mut out: Vec<(Vec<Base>, Vec<Variant>)> = Vec::new();
        for (reference, v) in pairs {
            match out.iter_mut().find(|(r, _)| r == &reference) {
                Some((_, vs)) => vs.push(v),
                None => out.push((reference, vec![v])),
            }
        }
        out
    }

    /// The sequence namespace separates otherwise-identical variants: the same
    /// edit on two assemblies must not collide.
    #[test]
    fn sequence_id_namespaces_the_id() {
        let reference = vec![base(0), base(1), base(0)];
        let v = Variant {
            pos: Interbase::new(1),
            del: vec![base(1)],
            ins: vec![base(0)],
        };
        let on_hg19 = variant_id("GRCh37:chr1", &reference, &v).unwrap();
        let on_hg38 = variant_id("GRCh38:chr1", &reference, &v).unwrap();
        assert_ne!(on_hg19, on_hg38);
    }

    /// Ids are stable and carry the documented prefix.
    #[test]
    fn id_is_stable_and_prefixed() {
        let reference = vec![base(0), base(0), base(1)];
        let v = Variant {
            pos: Interbase::new(0),
            del: vec![],
            ins: vec![base(0)],
        };
        let a = variant_id("seq", &reference, &v).unwrap();
        let b = variant_id("seq", &reference, &v).unwrap();
        assert_eq!(a, b);
        assert!(a.as_str().starts_with("cistron:va."));
        assert_eq!(a.as_str().len(), "cistron:va.".len() + 32);
    }

    /// A left-shiftable insertion and its already-canonical twin share an id.
    #[test]
    fn different_spellings_same_id() {
        // ref = A A C ; inserting an A can be written at pos 0, 1, or 2 — all the
        // same edit (an extra A in the A-run). Their ids must coincide.
        let reference = vec![base(0), base(0), base(1)];
        let at0 = Variant {
            pos: Interbase::new(0),
            del: vec![],
            ins: vec![base(0)],
        };
        let at1 = Variant {
            pos: Interbase::new(1),
            del: vec![],
            ins: vec![base(0)],
        };
        let at2 = Variant {
            pos: Interbase::new(2),
            del: vec![],
            ins: vec![base(0)],
        };
        let id0 = variant_id("seq", &reference, &at0).unwrap();
        let id1 = variant_id("seq", &reference, &at1).unwrap();
        let id2 = variant_id("seq", &reference, &at2).unwrap();
        assert_eq!(id0, id1);
        assert_eq!(id1, id2);
    }

    /// Pin the exact encoder output: id-equality tests pass under any consistent
    /// encoding, so the digest itself needs anchored vectors. RFC 4648 base64,
    /// with the URL-safe alphabet (`-`/`_`, never `+`/`/`).
    #[test]
    fn base64url_matches_known_vectors() {
        assert_eq!(base64url_nopad(b""), "");
        assert_eq!(base64url_nopad(b"f"), "Zg");
        assert_eq!(base64url_nopad(b"fo"), "Zm8");
        assert_eq!(base64url_nopad(b"foo"), "Zm9v");
        assert_eq!(base64url_nopad(b"foob"), "Zm9vYg");
        assert_eq!(base64url_nopad(b"fooba"), "Zm9vYmE");
        assert_eq!(base64url_nopad(b"foobar"), "Zm9vYmFy");
        assert_eq!(base64url_nopad(&[0xfa]), "-g"); // index 62 -> '-'
        assert_eq!(base64url_nopad(&[0xfe]), "_g"); // index 63 -> '_'
        assert_eq!(base64url_nopad(&[0xff, 0xff, 0xff]), "____");
    }

    #[test]
    fn display_equals_as_str() {
        let reference = vec![base(0), base(1), base(0)];
        let v = Variant {
            pos: Interbase::new(1),
            del: vec![base(1)],
            ins: vec![base(0)],
        };
        let id = variant_id("seq", &reference, &v).unwrap();
        assert_eq!(format!("{id}"), id.as_str());
        assert!(id.as_str().starts_with("cistron:va."));
    }

    #[test]
    fn ga4gh_bridge_produces_real_vrs_ids() {
        let acc = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";
        let reference = vec![base(0), base(1), base(2), base(3)]; // A C G T
                                                                  // SNV C>T at interbase 1 (already canonical): start=1, end=2, alt="T".
        let v = Variant {
            pos: Interbase::new(1),
            del: vec![base(1)],
            ins: vec![base(3)],
        };
        let id = ga4gh_allele_id(acc, &reference, &v).unwrap();
        assert_eq!(id, vrs::allele_id(acc, 1, 2, "T"));
        assert!(id.starts_with("ga4gh:VA."));
    }
}
