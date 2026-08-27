#!/bin/bash
# Task B3 production driver: chart-I CF-31 integer-slice sweep on aws-spot-11.
# Usage:   ./claude_z31_cf31_launch.sh [H] [LO] [HI]
# Default: H=30, shards [0,61) -- one shard per c4-index, sequential, resumable:
# a shard whose log already contains SWEEP_DONE is skipped, so after a spot
# interruption simply rerun this script on the same box (record: aws-spot-11).
# Sieve primes 1009 1013 1019 1021 (ascending -> empty-fiber early exit tests
# the cheapest prime first), CRT modulus M = 1.063e12, rational reconstruction
# bound 20000 (spurious rate ~1e-8/tuple), check prime 4001.
# DO NOT run this without the go-ahead: it is the production campaign.
set -u
cd "$(dirname "$0")"
H=${1:-30}
LO=${2:-0}
HI=${3:-$((2 * H + 1))}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-192}
BIN=./claude_z31_cf31_scan
for ((i = LO; i < HI; i++)); do
    LOG=$(printf "z31_cf31_int%02d_shard%02d.log" "$H" "$i")
    if grep -q SWEEP_DONE "$LOG" 2>/dev/null; then
        echo "skip shard $i (done)"
        continue
    fi
    echo "shard $i -> $LOG  ($(date -u +%FT%TZ) on $(hostname))"
    $BIN sweep "$H" 4 1009 1013 1019 1021 "$i" $((i + 1)) 20000 4001 int > "$LOG" 2>&1
done
echo "ALLSHARDS_DONE H=$H [$LO,$HI) $(date -u +%FT%TZ)"
