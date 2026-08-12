//! Domain enumeration and a reference `normalize`, plus the glue that hands
//! each `(raw -> out)` claim to the tlatools-rs oracle as a two-state walk.
//!
//! The real `cistron` crate does not exist yet. Until it does, `normalize`
//! here is a definitional brute force — the *leftmost parsimonious* variant
//! with the same denotation — which is correct by construction and makes the
//! wiring green. Swapping in an efficient `cistron::normalize` later means
//! deleting `normalize` below and re-pointing the test; the oracle keeps it
//! honest across the swap.

use serde_json::{json, Value as Json};

/// Bases are integers. Two letters is enough to exhibit repeats (hence
/// left-alignment); the spec itself is base-agnostic.
pub const ALPHABET: [i64; 2] = [0, 1];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    /// Interbase start: number of reference bases before the affected region.
    pub pos: usize,
    pub del: Vec<i64>,
    pub ins: Vec<i64>,
}

/// The alternate haplotype `v` produces against `reference`.
pub fn denotation(reference: &[i64], v: &Variant) -> Vec<i64> {
    let end = v.pos + v.del.len();
    let mut out = Vec::with_capacity(reference.len());
    out.extend_from_slice(&reference[..v.pos]);
    out.extend_from_slice(&v.ins);
    out.extend_from_slice(&reference[end..]);
    out
}

/// Well-formed iff in-bounds and `del` is exactly the reference at the locus.
pub fn well_formed(reference: &[i64], v: &Variant) -> bool {
    let end = v.pos + v.del.len();
    end <= reference.len() && v.del == reference[v.pos..end]
}

/// Blunt/parsimonious: del and ins share no leading and no trailing base.
pub fn parsimonious(v: &Variant) -> bool {
    if v.del.is_empty() || v.ins.is_empty() {
        return true;
    }
    v.del.first() != v.ins.first() && v.del.last() != v.ins.last()
}

/// Reference normalization: the leftmost parsimonious representative of the
/// edit `raw` denotes. Independent of the spec's local left-shift predicate,
/// so oracle agreement cross-validates that predicate rather than assuming it.
pub fn normalize(reference: &[i64], raw: &Variant) -> Variant {
    let target = denotation(reference, raw);
    let ins_max = reference.len() + raw.ins.len() + 1;

    let mut best: Option<Variant> = None;
    for pos in 0..=reference.len() {
        for dlen in 0..=(reference.len() - pos) {
            let del = reference[pos..pos + dlen].to_vec();
            for ins in all_seqs(ins_max) {
                let v = Variant {
                    pos,
                    del: del.clone(),
                    ins,
                };
                if denotation(reference, &v) == target && parsimonious(&v) {
                    // Order by (pos, dlen, ins): leftmost, then shortest deletion.
                    let better = match &best {
                        None => true,
                        Some(b) => (v.pos, v.del.len(), &v.ins) < (b.pos, b.del.len(), &b.ins),
                    };
                    if better {
                        best = Some(v);
                    }
                }
            }
        }
    }
    best.expect("every edit has a canonical representative")
}

/// A deliberately broken normalizer: returns the input untouched. Required to
/// *fail* the oracle on any non-canonical input, or the harness is vacuous.
pub fn normalize_broken(_reference: &[i64], raw: &Variant) -> Variant {
    raw.clone()
}

/// All base sequences of length `0..=max_len` over `ALPHABET`.
pub fn all_seqs(max_len: usize) -> Vec<Vec<i64>> {
    let mut out = vec![vec![]];
    let mut frontier = vec![vec![]];
    for _ in 0..max_len {
        let mut next = Vec::new();
        for seq in &frontier {
            for &b in &ALPHABET {
                let mut s = seq.clone();
                s.push(b);
                next.push(s);
            }
        }
        out.extend(next.iter().cloned());
        frontier = next;
    }
    out
}

/// Every reference of length `1..=max_ref_len` over `ALPHABET`.
pub fn all_refs(max_ref_len: usize) -> Vec<Vec<i64>> {
    all_seqs(max_ref_len)
        .into_iter()
        .filter(|r| !r.is_empty())
        .collect()
}

/// Every well-formed input variant on `reference`, deletions and insertions
/// bounded by `max_allele`.
pub fn all_variants(reference: &[i64], max_allele: usize) -> Vec<Variant> {
    let mut out = Vec::new();
    for pos in 0..=reference.len() {
        let dmax = (reference.len() - pos).min(max_allele);
        for dlen in 0..=dmax {
            let del = reference[pos..pos + dlen].to_vec();
            for ins in all_seqs(max_allele) {
                out.push(Variant {
                    pos,
                    del: del.clone(),
                    ins,
                });
            }
        }
    }
    out
}

fn variant_json(v: &Variant) -> Json {
    json!({ "pos": v.pos, "del": v.del, "ins": v.ins })
}

/// Build the oracle job for one `(raw -> out)` claim as a two-state walk.
pub fn job_json(spec: &str, reference: &[i64], raw: &Variant, out: &Variant) -> Json {
    let variant_schema = json!({
        "kind": "rec",
        "fields": {
            "pos": {"kind": "int"},
            "del": {"kind": "seq", "of": {"kind": "int"}},
            "ins": {"kind": "seq", "of": {"kind": "int"}}
        }
    });
    let root = json!({
        "ref": reference, "raw": variant_json(raw), "out": variant_json(raw), "done": false
    });
    let stepped = json!({
        "ref": reference, "raw": variant_json(raw), "out": variant_json(out), "done": true
    });
    json!({
        "spec": spec,
        "schema": {
            "ref": {"kind": "seq", "of": {"kind": "int"}},
            "raw": variant_schema,
            "out": variant_schema,
            "done": {"kind": "bool"}
        },
        "states": [root, stepped],
        "edges": [[0, 1, "normalize"]]
    })
}
