# cistron-identity

Content-addressed variant identifiers over [`cistron`](https://crates.io/crates/cistron)'s
normal form: two spellings of the same edit hash to the same id, so equality is
a hash compare (VRS-shaped truncated-SHA-512 digest).

Part of the [cistron workspace](https://github.com/copyleftdev/cistron). MIT OR Apache-2.0.
