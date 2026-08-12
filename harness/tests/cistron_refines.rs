//! The production algorithm under the same judge. `cistron::normalize` must (a)
//! refine `Normalize.tla` over the whole domain, and (b) agree base-for-base
//! with the harness's independent brute-force reference. (a) proves it obeys
//! the spec; (b) proves the spec's local left-shift rule equals the global
//! leftmost representative.

use cistron::{Base, Interbase, Variant as KVariant};
use cistron_oracle::{all_refs, all_variants, job_json, normalize as brute, Variant as HVariant};
use tla_oracle::{check, Job, Status};

const SPEC: &str = include_str!("../../specs/Normalize.tla");
const MAX_REF_LEN: usize = 4;
const MAX_ALLELE: usize = 2;

fn to_bases(xs: &[i64]) -> Vec<Base> {
    xs.iter()
        .map(|&x| Base::from_index(x as u8).expect("valid base index"))
        .collect()
}

fn from_bases(xs: &[Base]) -> Vec<i64> {
    xs.iter().map(|b| b.index() as i64).collect()
}

/// Run `cistron::normalize` and return the result in the harness's encoding.
fn cistron_normalize(reference: &[i64], raw: &HVariant) -> HVariant {
    let kv = KVariant {
        pos: Interbase::new(raw.pos),
        del: to_bases(&raw.del),
        ins: to_bases(&raw.ins),
    };
    let out = kv
        .normalize(&to_bases(reference))
        .expect("well-formed input normalizes");
    HVariant {
        pos: out.pos.get(),
        del: from_bases(&out.del),
        ins: from_bases(&out.ins),
    }
}

#[test]
fn cistron_normalize_refines_the_spec_and_matches_brute_force() {
    let mut checked = 0usize;
    for reference in all_refs(MAX_REF_LEN) {
        for raw in all_variants(&reference, MAX_ALLELE) {
            let out = cistron_normalize(&reference, &raw);

            // (b) agreement with the independent brute-force reference.
            let reference_answer = brute(&reference, &raw);
            assert_eq!(
                out, reference_answer,
                "cistron disagrees with brute force on ref={reference:?} raw={raw:?}"
            );

            // (a) the spec certifies the step.
            let job: Job = serde_json::from_value(job_json(SPEC, &reference, &raw, &out))
                .expect("job deserializes");
            let report = check(&job);
            assert_eq!(
                report.status,
                Status::Pass,
                "cistron output rejected: ref={reference:?} raw={raw:?} out={out:?}\n  {}",
                report.detail
            );
            checked += 1;
        }
    }
    eprintln!("cistron::normalize refined the spec on {checked} inputs");
    assert!(checked > 100, "domain collapsed to {checked}");
}
