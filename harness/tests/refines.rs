//! Enumerate a bounded genome domain, run `normalize` on every input, and make
//! the TLA+ spec certify each result. Coverage lives here, not in the oracle:
//! the oracle judges the walks we give it, so the numbers below *are* the claim.

use cistron_oracle::*;
use tla_oracle::{check, Job, Status};

const SPEC: &str = include_str!("../../specs/Normalize.tla");
const MAX_REF_LEN: usize = 4;
const MAX_ALLELE: usize = 2;

fn status(reference: &[i64], raw: &Variant, out: &Variant) -> tla_oracle::Report {
    let job: Job =
        serde_json::from_value(job_json(SPEC, reference, raw, out)).expect("job deserializes");
    check(&job)
}

/// Every `normalize` result over the domain must refine `Normalize.tla`.
#[test]
fn normalize_refines_the_spec() {
    let mut inputs = 0usize;
    let mut moved = 0usize; // inputs whose normal form differs from the input

    for reference in all_refs(MAX_REF_LEN) {
        for raw in all_variants(&reference, MAX_ALLELE) {
            let out = normalize(&reference, &raw);
            if out != raw {
                moved += 1;
            }
            let report = status(&reference, &raw, &out);
            assert_eq!(
                report.status,
                Status::Pass,
                "ref={reference:?} raw={raw:?} out={out:?}\n  {}",
                report.detail
            );
            inputs += 1;
        }
    }

    eprintln!("checked {inputs} inputs across the domain ({moved} were non-canonical)");
    assert!(inputs > 100, "domain collapsed to {inputs} inputs");
}

/// The canary: a normalizer that returns its input untouched must be *rejected*
/// on non-canonical inputs, and rejected at the edge (`Refines`), not the root.
#[test]
fn broken_normalizer_is_caught() {
    let mut rejected = 0usize;

    for reference in all_refs(MAX_REF_LEN) {
        for raw in all_variants(&reference, MAX_ALLELE) {
            let honest = normalize(&reference, &raw);
            if honest == raw {
                continue; // already canonical: identity is a correct answer here
            }
            let broken = normalize_broken(&reference, &raw);
            let report = status(&reference, &raw, &broken);
            assert_eq!(
                report.status,
                Status::Refines,
                "broken normalizer slipped past on ref={reference:?} raw={raw:?}\n  {}",
                report.detail
            );
            rejected += 1;
        }
    }

    eprintln!("canary rejected {rejected} non-canonical claims");
    assert!(rejected > 0, "canary never fired — the oracle is vacuous");
}
