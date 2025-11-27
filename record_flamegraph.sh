#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <output-prefix> <command> [args...]" >&2
  echo "Example: $0 autoc ~/GP/autoc/build/autoc -i ~/GP/autoc/autoc.ini" >&2
  exit 1
fi

OUTPUT_PREFIX="$1"
shift

PERF_BIN=${PERF_BIN:-perf}
PERF_FREQ=${PERF_FREQ:-997}
FLAMEGRAPH_DIR=${FLAMEGRAPH_DIR:-$HOME/FlameGraph}
STACK_COLLAPSE=${STACK_COLLAPSE:-$FLAMEGRAPH_DIR/stackcollapse-perf.pl}
FLAMEGRAPH_BIN=${FLAMEGRAPH_BIN:-$FLAMEGRAPH_DIR/flamegraph.pl}

if ! command -v "$PERF_BIN" >/dev/null 2>&1; then
  echo "Error: $PERF_BIN not found in PATH" >&2
  exit 1
fi

if [[ ! -x "$STACK_COLLAPSE" || ! -x "$FLAMEGRAPH_BIN" ]]; then
  cat >&2 <<EOF
Error: FlameGraph utilities not found.
Set FLAMEGRAPH_DIR to the directory containing stackcollapse-perf.pl and flamegraph.pl.
Repo: https://github.com/brendangregg/FlameGraph
EOF
  exit 1
fi

PERF_DATA="${OUTPUT_PREFIX}.perf.data"
FOLDED="${OUTPUT_PREFIX}.perf.folded"
SVG="${OUTPUT_PREFIX}.perf.svg"

echo "[perf] recording to ${PERF_DATA}"
sudo --preserve-env=PERF_BIN "$PERF_BIN" record \
  --call-graph dwarf \
  -F "$PERF_FREQ" \
  -o "$PERF_DATA" \
  -- "$@"

echo "[perf] converting to folded stack ${FOLDED}"
"$PERF_BIN" script -i "$PERF_DATA" | "$STACK_COLLAPSE" > "$FOLDED"

echo "[FlameGraph] generating ${SVG}"
"$FLAMEGRAPH_BIN" "$FOLDED" > "$SVG"

echo "Flame graph written to ${SVG}"
