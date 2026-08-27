#!/usr/bin/env bash
set -euo pipefail

# Standalone driver for the bounded exact a=1/e scan.  The Magma
# continuation needs the reconstruction variables but the reconstruction
# source normally terminates; this driver composes them without requiring
# the user to edit either source file.

height="${1:-200}"
prime_bound="${2:-43}"
timeout_seconds="${SCAN_TIMEOUT_SECONDS:-300}"

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
cd "$repo_dir"

sed -e '$r code/contact6_m612_relative3_rational_a_scan_continuation.m' \
    -e '$d' code/contact6_m612_relative3_exact_reconstruct.m \
| timeout "$timeout_seconds" magma -b do_primitive24:=false \
    print_maps:=false height:="$height" prime_bound:="$prime_bound"
