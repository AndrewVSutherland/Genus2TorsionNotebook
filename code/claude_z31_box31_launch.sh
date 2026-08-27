#!/bin/bash
# Task B5 production launcher for claude_z31_box31 on a 192-core spot box.
# NOT run automatically -- launch by hand after checking core availability
# (spot-11 is shared: ps aux --sort=-%cpu | head first, coordinate with any
# sibling scans, and reduce NSH*THR to the free-core budget if needed).
#
# Usage (on the spot box, from ~/z31):
#   ./claude_z31_box31_launch.sh [H] [NSHARDS] [THREADS_PER_SHARD]
# Defaults: H=43, 16 shards x 12 threads = 192 cores, ~21.5 h wall at the
# measured 1.253M curves/s/core (1.86e13 tested models after dedupe).
#
# Each shard writes results/box31_H{H}_s{S}of{N}.log with SURV / PROGRESS /
# BOX_DONE markers.  After a spot interruption, rerun ONLY the shards whose
# log lacks a BOX_DONE line (same command line; shards are independent and
# idempotent).  Survivors:  grep -h '^SURV' results/box31_H*_s*.log
# go to the exact stage:  magma -b C6:=.. C5:=.. C4:=.. C3:=.. C2:=.. C1:=.. C0:=.. \
#                           claude_z31_box31_survcheck.m
H=${1:-43}
NSH=${2:-16}
THR=${3:-12}
cd "$(dirname "$0")"
mkdir -p results
echo "launching box H=$H in $NSH shards x $THR threads on $(hostname)"
for s in $(seq 0 $((NSH - 1))); do
    log=results/box31_H${H}_s${s}of${NSH}.log
    if grep -q BOX_DONE "$log" 2>/dev/null; then
        echo "shard $s already done, skipping"
        continue
    fi
    if [ -f "$log" ] && [ -n "$(find "$log" -mmin -3 2>/dev/null)" ]; then
        echo "shard $s log active in the last 3 min -- probably still RUNNING, skipping (rm log to force)"
        continue
    fi
    OMP_NUM_THREADS=$THR nohup ./claude_z31_box31 box "$H" "$s" "$NSH" > "$log" 2>&1 &
    echo "shard $s pid $!"
done
echo "monitor:  grep -h PROGRESS results/box31_H${H}_s*.log | tail; grep -hc SURV results/box31_H${H}_s*.log"
