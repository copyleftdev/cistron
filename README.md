# cistron

[![CI](https://github.com/copyleftdev/cistron/actions/workflows/ci.yml/badge.svg)](https://github.com/copyleftdev/cistron/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license)
[![tip](https://img.shields.io/badge/tip-%40copyleftdev-ff69b4.svg)](https://tokentip.to/@copyleftdev)

The core primitive for the genome and the rules that apply to it: a
convention-safe *located allele* plus the normalization and equivalence rules
that decide when two spellings are the same variant.

Normalization is validated **byte-for-byte against `bcftools norm`** over
291,600 cases (see the differential below); the kernel and every boundary crate
are mutation-tested to **0 survivors**, and both parsers are fuzzed.

Nothing here is ML or annotation. It is the atom every layer above imports and
the small set of rules that make that atom trustworthy — coordinate-correct,
denotation-preserving, and content-addressable.

## Layout

```
cistron/          the kernel crate: Base, Interbase, Variant, normalize
identity/         cistron-identity: content-addressed VariantId over the normal form
vcf/              cistron-vcf: 1-based anchored VCF records <-> interbase blunt variants
hgvs/             cistron-hgvs: 3'-shifted genomic (g.) nomenclature <-> interbase variants
specs/            the rules as TLA+, with an independent TLC cross-check
harness/          runs cistron::normalize against the spec via the tlatools-rs
                  oracle over an enumerated genome domain (2058 inputs)
```

## The idea

The rule is written once, in TLA+, and everything else is held to it:

- **`specs/Normalize.tla`** *recognises* a correct normal form (it does not
  recompute one). The oracle (github.com/copyleftdev/tlatools-rs) is a
  refinement judge: the harness runs `cistron::normalize`, hands each
  `raw → out` step in as a two-state walk, and the spec certifies it. Coverage
  lives in the harness, so the printed input count *is* the claim.
- **`specs/NormalizeCheck.tla`** re-checks the same semantics with **real TLC**,
  to guard against the oracle grading its own homework. It proves the two
  properties that make a variant content-addressable:
  - **Sound** — the canonical form of an edit is *unique* (same denotation +
    both canonical ⇒ identical).
  - **Total** — every variant *has* a canonical form.
- **`specs/Coords.tla`** checks the coordinate conventions: interbase ↔ 1-based
  round-trips, and half-open overlap agrees with inclusive-space overlap.

Every TLC spec ships a **canary** invariant that must be *violated* (e.g.
left-alignment dropped, or half-open overlap written with `<=`), so a green run
proves the check is not vacuous.

## Convention

LEFT-aligned + blunt/parsimonious — the bcftools-norm / GA4GH-VRS convention.
Coordinates are **interbase** (`pos` = reference bases before the locus,
0-based); the 0/1-based edge conventions are converted at the boundary via
`Interbase::from_zero_based_start` / `from_one_based`. HGVS flips the shift to
the RIGHT and is a separate, deliberately absent module.

## Rules enforced

- **Reference-agreement** — `del` must equal the reference at its locus
  (`Variant::check`); the largest class of real variant bugs.
- **Preservation** — `normalize` never changes what a variant denotes.
- **Canonicity** — the result is parsimonious and cannot shift left.
- **Idempotence** — falls out for free (re-normalizing is a stutter step).
- **Uniqueness / totality** — proven by TLC, above.

## Running it

```
cargo test              # kernel property tests + oracle refinement over 2058 inputs
specs/check.sh          # independent TLC cross-check (real invariants + canaries)
```

Rigor beyond the internal validators:

```
cargo mutants -p cistron            # kernel mutation score: 0 missed
cargo +nightly fuzz run hgvs_parse  # parser panic-safety (found+fixed an overflow)
cargo +nightly fuzz run vcf_parse
BCFTOOLS=/path/to/bcftools cargo run -p cistron-difftest
    # EXTERNAL oracle: cistron::normalize vs `bcftools norm`, byte-for-byte.
    # 291,600 records agree (79,632 required realignment) — the convention is
    # validated against the tool the community trusts, not just our own spec.
```

Representative output:

```
cistron::normalize refined the spec on 714 inputs
canary rejected 376 non-canonical claims
all cross-checks passed
```

## Status and next

Phase 0 is complete: the located-allele type, a validated left-align
`normalize`, and the spec + cross-check that pin it. Bounds are `MAX_REF_LEN=3`,
`MAX_ALLELE=2` over a 2-letter alphabet — enough to exhibit repeats and
left-alignment; widen in the harness tests to scale the exhaustive pass.

**`identity`** is built: `variant_id(sequence_id, reference, variant)` normalizes
then digests `(sequence_id, start, end, alt)` — VRS's algorithm (truncated
SHA-512, 24 bytes, base64url), `cistron:va.` prefix. The identity law
(`id equal ⟺ same denotation`, on one sequence) is exhaustively tested, so
equality is a hash compare. Byte-for-byte GA4GH-VRS compat is a documented swap
of the serialization only.

**`cistron-vcf`** is built: `to_variant`/`parse_line` map VCF's 1-based anchored
records into interbase blunt variants (the anchor base trims away under
`normalize`); `to_vcf` emits, re-introducing an anchor from the reference
(left-anchored, right-anchored only at the reference start). Round-trip is
proven lossless up to normalization — emit-then-parse lands on the same normal
form — and no-op / whole-reference-start-deletion cases are reported, not
dropped. Symbolic ALTs, breakends, `*`, and non-ACGT bases are explicit errors.

**`cistron-hgvs`** is built. It exercised the design's hardest claim — a second,
*opposite* convention — and it held: the kernel gained `normalize_right` (the
3'-most dual of `normalize`), and HGVS renders/parses genomic `g.` expressions
(`>`, `del`, `ins`, `dup`, `delins`), including the tandem-`dup`-not-`ins` rule
and the fact that HGVS omits deleted/duplicated bases (so parsing needs the
reference). Round-trip render→parse recovers the left canonical over the domain.
A key thing surfaced: in a tandem repeat the inserted sequence *rotates* between
the left and right forms (`AC`↔`CA`) — same denotation, different string —
which is exactly *why* VCF and HGVS disagree, and the kernel now models it.

Next layers, each a crate *depending on* the kernel, never the reverse:

- **liftover** — the one remaining coordinate rule: fallible, never silent.
- **VRS byte-compat** — swap `identity::canonical_bytes` for VRS `ga4gh_serialize`
  to emit real `ga4gh:VA.` ids.

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at
your option. Unless you explicitly state otherwise, any contribution
intentionally submitted for inclusion in the work by you shall be dual licensed
as above, without any additional terms or conditions.
