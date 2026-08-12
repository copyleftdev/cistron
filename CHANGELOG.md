# Changelog

All notable changes are documented here. This project follows
[Semantic Versioning](https://semver.org). Releases are cut by
[release-plz](https://release-plz.dev) from conventional-commit history; use
`scripts/bump.sh` for a manual bump.

## [Unreleased]

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
