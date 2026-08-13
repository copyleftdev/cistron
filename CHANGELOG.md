# Changelog

All notable changes are documented here. This project follows
[Semantic Versioning](https://semver.org). Releases are cut by
[release-plz](https://release-plz.dev) from conventional-commit history; use
`scripts/bump.sh` for a manual bump.

## [Unreleased]

- `cistron-hgvs` — genomic rendering validated byte-for-byte against biocommons
  `hgvs` (2,000-vector differential); added inversion (`inv`) detection/parsing
  and start-of-sequence `delins`, the two rules the differential surfaced.

- `cistron-identity::vrs` — GA4GH `ReferenceLengthExpression` support; `ga4gh_allele_id`
  now matches vrs-python byte-for-byte for large-repeat variants too (literal +
  RLE state decision), validated by a 2,000-vector differential.

- `Variant::fully_justified` — GA4GH-VRS fully-justified normalization, matching
  vrs-python byte-for-byte (2,000-vector differential); `ga4gh_allele_id` now
  uses it, so VRS ids match the reference for tandem-repeat indels too.

- `cistron-identity::vrs` — real GA4GH VRS identifiers (`ga4gh:VA.`/`SL.`/`SQ.`),
  byte-for-byte with vrs-python (validated against the spec example and a
  committed 2,001-vector reference corpus).

- `cistron-liftover` — fallible, never-silent interbase coordinate liftover
  (chain algebra + UCSC `.chain` file parser; `Unmapped`/`Split` errors,
  strand-aware, fuzzed and mutation-clean).

## [0.1.0]

Initial release.

- `cistron` — the located-allele kernel: interbase coordinates, reference
  agreement, and left-aligned parsimonious `normalize` (plus `normalize_right`).
  Validated byte-for-byte against `bcftools norm` over 291,600 records and
  mutation-tested to zero surviving mutants.
- `cistron-identity` — content-addressed variant identifiers over the normal
  form (VRS-shaped digest).
- `cistron-vcf` — the VCF boundary (1-based anchored records ↔ interbase).
- `cistron-hgvs` — the HGVS genomic (`g.`) boundary (3'-shifted nomenclature).
