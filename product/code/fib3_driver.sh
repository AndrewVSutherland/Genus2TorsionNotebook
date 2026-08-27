#!/bin/bash
# fib3_driver.sh — crash-resuming driver for lane_2266_sigma_fib2.m
# (Magma segfaults state-dependently after many descents; a segfault is not
# catchable in-language, so drive from outside: on abnormal exit, retry the
# crashing fiber once in a fresh process; if it crashes twice in a row, mark
# it CRASHSKIP and move past it.  Totals are recomputed from the full log.)
# Usage: ./fib3_driver.sh <Sig> [Skip0]
set -u
SIG=$1
SKIP=${2:-0}
LOG=../logs/sigfib3_si${SIG}.log
PROG=../logs/sigmaF${SIG}b.progress
PREVCRASH=-1
for attempt in $(seq 1 60); do
  timeout --signal=KILL 9000 magma -b Sig:=$SIG Skip:=$SKIP lane_2266_sigma_fib2.m >> $LOG 2>&1
  rc=$?
  if tail -40 "$LOG" | grep -q "SEARCH_DONE 2266sigfib2"; then
    echo "DRIVER: natural end (attempt $attempt, rc=$rc)" >> $LOG
    exit 0
  fi
  last=$(tail -1 "$PROG" | sed -n 's/.* FIB \([0-9]*\)\/.*/\1/p')
  if [ -z "$last" ]; then echo "DRIVER: no progress line, aborting" >> $LOG; exit 1; fi
  if [ "$last" = "$PREVCRASH" ]; then
    echo "CRASHSKIP fiber $last (Sig=$SIG, crashed twice)" >> $LOG
    SKIP=$last
    PREVCRASH=-1
  else
    echo "DRIVER: crash at fiber $last (rc=$rc), retrying once fresh" >> $LOG
    SKIP=$((last-1))
    PREVCRASH=$last
  fi
done
echo "DRIVER: attempt budget exhausted" >> $LOG
exit 1
