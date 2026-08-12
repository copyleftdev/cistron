# cistron-identity

Content-addressed variant identifiers over [`cistron`](https://crates.io/crates/cistron)'s
normal form.

- `variant_id` — a compact internal content-address (`cistron:va.`); equality is
  a hash compare (identity law exhaustively tested).
- `vrs` — **real GA4GH VRS identifiers** (`ga4gh:VA.`, `ga4gh:SL.`, `ga4gh:SQ.`),
  **byte-for-byte compatible with vrs-python** (validated against the VRS 2.0
  spec example and a committed 2,001-vector reference corpus).

Part of the [cistron workspace](https://github.com/copyleftdev/cistron). MIT OR Apache-2.0.
