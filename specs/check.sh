#!/usr/bin/env bash
# Independent TLC cross-check of the spec semantics. Each spec ships a real
# invariant that must HOLD and a canary invariant that must be VIOLATED, so a
# green run proves the check is non-vacuous. Exits nonzero on any surprise.
set -euo pipefail
cd "$(dirname "$0")"

JAR="${TLA2TOOLS:-/home/ops/Project/tla-for-ai/vendor/tla2tools.jar}"
# TLC exits nonzero when it finds a violation, so tolerate that and inspect the
# text — the whole point of the canary is a run that "fails".
tlc() { java -cp "$JAR" tlc2.TLC -deadlock -cleanup -metadir "$(mktemp -d)" -config "$1" "$2" 2>&1 || true; }

expect_holds() { # cfg module
  if tlc "$1" "$2" | grep -q "No error has been found"; then
    echo "  ok: $1 holds"
  else
    echo "  FAIL: $1 was expected to hold"; exit 1
  fi
}
expect_violated() { # cfg module
  if tlc "$1" "$2" | grep -qiE "is violated|equal to FALSE"; then
    echo "  ok: $1 violated (canary fired)"
  else
    echo "  FAIL: $1 canary did not fire"; exit 1
  fi
}

echo "Normalize semantics (Sound /\\ Total):"
expect_holds    NormalizeCheck.cfg NormalizeCheck.tla
expect_violated Canary.cfg         NormalizeCheck.tla

echo "Coordinate conventions (round-trip + overlap agreement):"
expect_holds    Coords.cfg         Coords.tla
expect_violated CoordsCanary.cfg   Coords.tla

echo "all cross-checks passed"
