//! Golden vectors pinning each branch of the VRS state decision (literal vs
//! run-length expression) with exact vrs-python identifiers, so the branch logic
//! and the repeat-subunit computation are anchored, not just covered in
//! aggregate by the random differential.

use cistron::{Base, Interbase, Variant};
use cistron_identity::ga4gh_allele_id;

const ACC: &str = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl";

fn bases(s: &str) -> Vec<Base> {
    s.chars()
        .map(|c| match c {
            'A' => Base::A,
            'C' => Base::C,
            'G' => Base::G,
            'T' => Base::T,
            other => panic!("bad base {other}"),
        })
        .collect()
}

fn id(seq: &str, start: usize, end: usize, alt: &str) -> String {
    let reference = bases(seq);
    let raw = Variant {
        pos: Interbase::new(start),
        del: reference[start..end].to_vec(),
        ins: bases(alt),
    };
    ga4gh_allele_id(ACC, &reference, &raw).unwrap()
}

#[test]
fn vrs_state_decision_golden_ids() {
    // substitution / complex edit -> LiteralSequenceExpression
    assert_eq!(
        id("GAAAT", 1, 4, "CC"),
        "ga4gh:VA.Sq6v2FRiY6CTyZ8wUHdNfFxRYSsF4xps"
    );
    // non-repeat insertion (reference not expanded) -> literal
    assert_eq!(
        id("ACGT", 1, 1, "TT"),
        "ga4gh:VA.pZhTBnmzBxKdz6n7PWts17U9hmpUWwQY"
    );
    // non-repeat deletion -> RLE with length 0, repeatSubunitLength 1
    assert_eq!(
        id("ACGT", 1, 2, ""),
        "ga4gh:VA.8STe303q5gpHLvwcTJUq5oXBdewPEAmv"
    );
    // dinucleotide-repeat insertion -> RLE, repeatSubunitLength = seed (2)
    assert_eq!(
        id("GATATATATATATATATC", 1, 1, "AT"),
        "ga4gh:VA.CWU-PtkMz4ija29luEZSa-opbY5yz41s"
    );
    // insertion whose subunit is a PROPER divisor of the seed (seed 4 -> subunit 1):
    // exercises factors_desc reaching a later factor.
    assert_eq!(
        id("ATGTCGTAGGCCC", 4, 4, "TTTT"),
        "ga4gh:VA.LMmbQaZ1LLNrvEY-SrtFojDPYvUjOYvI"
    );
    // insertion where the subunit is a proper factor > 1 (seed 4 -> subunit 2):
    // factors_desc's factor `2` must be present and selected, or the id changes.
    assert_eq!(
        id("AAACCAA", 1, 1, "AAAA"),
        "ga4gh:VA.EnMi_YGcQNgr4aODGOtJBmLC_8z3ue0U"
    );
}
