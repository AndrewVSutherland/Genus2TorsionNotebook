#!/bin/bash
# box_driver.sh — crash/TimeCap-resuming driver for the CORRECTED box fleet
# (lane_2266_sigma_fib2.m with verified saturation-adoption), logging to
# sigfib4_si<k>.log.  Usage: ./box_driver.sh <Sig> [Skip0]
set -u
SIG=$1
SKIP=${2:-0}
LOG=../logs/sigfib4_si${SIG}.log
PROG=../logs/sigmaF${SIG}b.progress
PREVCRASH=-1
# fresh start (Skip0=0): truncate the log so stale dispositions from earlier
# fleets cannot bleed into accounting or regenerated inputs (codex); resumes
# (Skip0>0) append.
if [ "$SKIP" = "0" ]; then
  : > $LOG
fi
for attempt in $(seq 1 80); do
  timeout --signal=KILL 12000 magma -b Sig:=$SIG Skip:=$SKIP TimeCap:=10800 lane_2266_sigma_fib2.m >> $LOG 2>&1
  rc=$?
  seg_done=$(tail -50 "$LOG" | grep -c "SEARCH_DONE 2266sigfib2")
  seg_cap=$(tail -50 "$LOG" | grep -oE "TIMECAP [0-9]+ s at fiber [0-9]+" | tail -1 | grep -oE "fiber [0-9]+" | grep -oE "[0-9]+")
  if [ "$seg_done" -ge 1 ] && [ -z "$seg_cap" ]; then
    echo "DRIVER: natural end (attempt $attempt, rc=$rc)" >> $LOG
    exit 0
  fi
  if [ -n "$seg_cap" ]; then
    echo "DRIVER: TimeCap at fiber $seg_cap, resuming" >> $LOG
    SKIP=$((seg_cap-1))
    PREVCRASH=-1
    continue
  fi
  last=$(tail -1 "$PROG" | sed -n 's/.* FIB \([0-9]*\)\/.*/\1/p')
  if [ -z "$last" ]; then echo "DRIVER: no progress line, aborting" >> $LOG; exit 1; fi
  if [ "$last" = "$PREVCRASH" ]; then
    echo "CRASHSKIP fiber $last (Sig=$SIG, crashed twice)" >> $LOG
    SKIP=$last
    PREVCRASH=-1
  else
    echo "DRIVER: crash/kill at fiber $last (rc=$rc), retrying once fresh" >> $LOG
    SKIP=$((last-1))
    PREVCRASH=$last
  fi
done
echo "DRIVER: attempt budget exhausted" >> $LOG
exit 1
