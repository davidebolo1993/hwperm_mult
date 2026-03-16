#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT_DIR/src/hwperm_mult.cpp"
OUT="$ROOT_DIR/hwperm_mult"

echo "Compiling $SRC -> $OUT"
g++ -O3 -std=c++17 -Wall -Wextra -pedantic "$SRC" -o "$OUT"
echo "Done."
