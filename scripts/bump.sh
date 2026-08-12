#!/usr/bin/env bash
# Bump the shared workspace version so a release can never collide with an
# already-published one. Every crate inherits `version.workspace`, and the
# internal dependency requirements live in `[workspace.dependencies]`, so all
# versions are in the root Cargo.toml and move together in one edit.
#
#   scripts/bump.sh          # patch: 0.1.0 -> 0.1.1
#   scripts/bump.sh minor    # 0.1.0 -> 0.2.0
#   scripts/bump.sh major    # 0.1.0 -> 1.0.0
set -euo pipefail
cd "$(dirname "$0")/.."

part="${1:-patch}"
cur=$(grep -m1 -E '^version = "[0-9]+\.[0-9]+\.[0-9]+"' Cargo.toml |
    sed -E 's/.*"([0-9]+\.[0-9]+\.[0-9]+)".*/\1/')
[ -n "$cur" ] || {
    echo "could not find workspace version in Cargo.toml" >&2
    exit 1
}

IFS=. read -r ma mi pa <<<"$cur"
case "$part" in
major) ma=$((ma + 1)); mi=0; pa=0 ;;
minor) mi=$((mi + 1)); pa=0 ;;
patch) pa=$((pa + 1)) ;;
*)
    echo "usage: $0 [major|minor|patch]" >&2
    exit 1
    ;;
esac
new="$ma.$mi.$pa"

# Every "$cur" in the root manifest is ours (workspace.package +
# workspace.dependencies); external deps and rust-version never match X.Y.Z.
sed -i "s/\"$cur\"/\"$new\"/g" Cargo.toml
echo "workspace version: $cur -> $new"
