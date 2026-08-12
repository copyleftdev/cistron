# cistron-liftover

Fallible, never-silent interbase coordinate liftover between assemblies. A
coordinate that cannot map cleanly returns a typed error (`Unmapped`/`Split`) —
never a plausible-but-wrong number. The chain *algebra*; UCSC chain-file parsing
is a boundary concern layered on top.

Part of the [cistron workspace](https://github.com/copyleftdev/cistron). MIT OR Apache-2.0.
