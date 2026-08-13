//! GA4GH VRS computed identifiers, byte-for-byte with the specification.
//!
//! These reproduce the identifiers other VRS tools produce, so `cistron` ids
//! interoperate. The algorithm is validated against the official worked example
//! (see the tests): a SequenceLocation digest of `wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz`
//! and an Allele id of `ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt`.
//!
//! Serialization follows the VRS rule set: keys ordered by Unicode code point,
//! compact separators, nested *identifiable* objects (the SequenceLocation
//! inside an Allele) replaced by their digest, and only GA4GH-Digest-Keys
//! fields included. The `SequenceReference` is a reference, not identifiable, so
//! it is serialized inline.

use crate::base64url_nopad;
use sha2::{Digest, Sha512};

/// The GA4GH truncated digest: base64url of the leftmost 24 bytes of SHA-512.
pub fn sha512t24u(blob: &[u8]) -> String {
    let digest = Sha512::digest(blob);
    base64url_nopad(&digest[..24])
}

/// The refget digest of a sequence (case-normalized to upper), as a bare digest.
pub fn sequence_digest(sequence: &[u8]) -> String {
    sha512t24u(&sequence.to_ascii_uppercase())
}

/// The refget sequence identifier: `ga4gh:SQ.<digest>`.
pub fn sequence_id(sequence: &[u8]) -> String {
    format!("ga4gh:SQ.{}", sequence_digest(sequence))
}

/// Canonical serialization of a SequenceLocation. `refget_accession` carries its
/// own `SQ.` prefix (a refget accession), not the `ga4gh:` CURIE prefix.
fn location_serialize(refget_accession: &str, start: u64, end: u64) -> String {
    format!(
        "{{\"end\":{end},\
         \"sequenceReference\":{{\"refgetAccession\":\"{refget_accession}\",\"type\":\"SequenceReference\"}},\
         \"start\":{start},\"type\":\"SequenceLocation\"}}"
    )
}

/// The bare digest of a SequenceLocation.
pub fn location_digest(refget_accession: &str, start: u64, end: u64) -> String {
    sha512t24u(location_serialize(refget_accession, start, end).as_bytes())
}

/// The SequenceLocation identifier: `ga4gh:SL.<digest>`.
pub fn location_id(refget_accession: &str, start: u64, end: u64) -> String {
    format!("ga4gh:SL.{}", location_digest(refget_accession, start, end))
}

/// Canonical serialization of an Allele: the location is replaced by its digest.
fn allele_serialize(location_digest: &str, alt: &str) -> String {
    format!(
        "{{\"location\":\"{location_digest}\",\
         \"state\":{{\"sequence\":\"{alt}\",\"type\":\"LiteralSequenceExpression\"}},\
         \"type\":\"Allele\"}}"
    )
}

/// The bare digest of an Allele at `[start, end)` on `refget_accession` with
/// alternate sequence `alt` (interbase coordinates; `alt` is IUPAC, uppercase).
pub fn allele_digest(refget_accession: &str, start: u64, end: u64, alt: &str) -> String {
    let loc = location_digest(refget_accession, start, end);
    sha512t24u(allele_serialize(&loc, alt).as_bytes())
}

/// The Allele identifier: `ga4gh:VA.<digest>`.
pub fn allele_id(refget_accession: &str, start: u64, end: u64, alt: &str) -> String {
    format!(
        "ga4gh:VA.{}",
        allele_digest(refget_accession, start, end, alt)
    )
}

/// Canonical serialization of an Allele whose state is a ReferenceLengthExpression
/// (a large repeat expansion). The `sequence` field is *excluded* from the digest,
/// so the id depends only on the interval, `length`, and `repeatSubunitLength`.
fn rle_allele_serialize(location_digest: &str, length: u64, repeat_subunit_length: u64) -> String {
    format!(
        "{{\"location\":\"{location_digest}\",\
         \"state\":{{\"length\":{length},\"repeatSubunitLength\":{repeat_subunit_length},\
         \"type\":\"ReferenceLengthExpression\"}},\"type\":\"Allele\"}}"
    )
}

/// The Allele identifier for a run-length-expressed variant (`ga4gh:VA.…`).
pub fn rle_allele_id(
    refget_accession: &str,
    start: u64,
    end: u64,
    length: u64,
    repeat_subunit_length: u64,
) -> String {
    let loc = location_digest(refget_accession, start, end);
    format!(
        "ga4gh:VA.{}",
        sha512t24u(rle_allele_serialize(&loc, length, repeat_subunit_length).as_bytes())
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const ACC: &str = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";

    /// The entire pipeline, from structured inputs, against the VRS 2.0 spec's
    /// worked example (NC_000019.10:g.44908822 A>T style SNV).
    #[test]
    fn official_example_end_to_end() {
        assert_eq!(
            location_digest(ACC, 44908821, 44908822),
            "wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz"
        );
        assert_eq!(
            location_id(ACC, 44908821, 44908822),
            "ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz"
        );
        assert_eq!(
            allele_digest(ACC, 44908821, 44908822, "T"),
            "0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
        );
        assert_eq!(
            allele_id(ACC, 44908821, 44908822, "T"),
            "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
        );
    }

    /// Run-length-expression Allele ids against vrs-python reference values
    /// (a 4-base deletion and a 10-base insertion in a 60-base A-run at [1,61)).
    #[test]
    fn rle_allele_ids() {
        assert_eq!(
            location_digest(ACC, 1, 61),
            "L3WZOuLPUmIUCJ-hKYObXR_LieL7oXSB"
        );
        assert_eq!(
            rle_allele_id(ACC, 1, 61, 56, 4),
            "ga4gh:VA.2QZjIh4_0ipE8mdtFd1RGOAM8wElTaFj"
        );
        assert_eq!(
            rle_allele_id(ACC, 1, 61, 70, 10),
            "ga4gh:VA.v9pAnGWFLKu-ca-gFBGXHJBsL7P5ZmgL"
        );
    }

    /// sha512t24u produces 32 unpadded base64url chars over the 24-byte digest.
    #[test]
    fn digest_shape() {
        let d = sha512t24u(b"");
        assert_eq!(d.len(), 32);
        assert!(!d.contains('=') && !d.contains('+') && !d.contains('/'));
    }

    /// Refget sequence digests against ga4gh.core.sha512t24u vectors.
    #[test]
    fn refget_sequence_digests() {
        assert_eq!(sequence_digest(b"ACGT"), "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2");
        assert_eq!(sequence_digest(b""), "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc");
        // residues are upper-cased before digesting.
        assert_eq!(sequence_digest(b"acgt"), sequence_digest(b"ACGT"));
        assert_eq!(
            sequence_id(b"ACGT"),
            "ga4gh:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
        );
    }

    /// A deletion has an empty alternate sequence; the serialization must still
    /// be well-formed.
    #[test]
    fn empty_alt_is_a_deletion() {
        let id = allele_id(ACC, 100, 105, "");
        assert!(id.starts_with("ga4gh:VA."));
        assert_eq!(id.len(), "ga4gh:VA.".len() + 32);
    }
}
