#!/bin/bash
# claude_magma_slot.sh -- run Magma while holding one of N GLOBAL slots.
#
# Why this exists: on 2026-07-25 nine autonomous lanes were each told "at most one
# Magma job, wait if the machine already has >= 12".  Nine agents each judging
# "one job is fine" does not compose into a global limit -- claude-box reached
# load 70 with 38 concurrent magma.exe on 32 vCPUs.  Instruction is not a
# mechanism.  This wrapper is the mechanism.
#
# Usage (drop-in replacement for `magma`):
#     code/claude_magma_slot.sh -b MYVAR:=1 code/myscript.m > results/my.log 2>&1 &
# Env:
#     MAGMA_SLOTS   number of concurrent slots (default 10)
#     SLOT_TIMEOUT  seconds to wait for a slot before giving up (default 86400)
#
# The lock is held on a file descriptor that is inherited by the exec'd magma
# process, so the slot is released automatically when Magma exits -- including
# on crash, kill, or an agent dying mid-run.  No cleanup required, no stale locks.

SLOTS=${MAGMA_SLOTS:-10}
LOCKDIR=${MAGMA_LOCKDIR:-/tmp/claude-magma-slots}
TIMEOUT=${SLOT_TIMEOUT:-86400}
mkdir -p "$LOCKDIR" 2>/dev/null

start=$(date +%s)
while :; do
  for i in $(seq 1 "$SLOTS"); do
    exec 9>"$LOCKDIR/slot$i" 2>/dev/null || continue
    if flock -n 9; then
      # fd 9 stays open across exec -> the slot is held for the life of Magma
      exec magma "$@"
    fi
    exec 9>&-
  done
  now=$(date +%s)
  if [ $(( now - start )) -ge "$TIMEOUT" ]; then
    echo "claude_magma_slot: no slot within ${TIMEOUT}s (SLOTS=$SLOTS)" >&2
    exit 75
  fi
  sleep 5
done
