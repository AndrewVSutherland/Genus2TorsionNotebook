#!/bin/bash
# (8,8) lane, step 1.4 wrapper: sequential per-base exact-J1 probes.
# One Magma process at a time (occupies ONE of the 3 allowed Magma slots).
# External timeout per base: MordellWeilGroup can hang on hard bases, and a
# killed -b run leaves an empty log -- the master log records TIMEOUT markers.
# Base (3,1/3) is skipped (already done: data/claude_prod_09_88_double_stage1.txt).
#
# Run: nohup ./code/claude_prod_09_88_j1exact_run.sh > /dev/null 2>&1 &

cd torsion_jac || exit 1
MASTER=results/claude_prod_09_88_j1exact_master.log

# m_num/m_den:n_num/n_den  -- (3,5) first (known rank 3), (2,1) last (known slow)
BASES="3/1:5/1 2/1:3/1 3/1:2/1 2/1:5/1 3/1:1/1 4/1:1/1 1/2:2/1 2/3:1/1 5/1:2/1 2/1:1/3 3/1:4/1 2/1:1/1"

echo "J1EXACT_SWEEP_START $(date '+%F %T')" >> "$MASTER"
for bp in $BASES; do
  mm=${bp%%:*}; nn=${bp##*:}
  mn=${mm%%/*}; md=${mm##*/}
  nnu=${nn%%/*}; nd=${nn##*/}
  log=results/claude_prod_09_88_j1exact_m${mn}_${md}_n${nnu}_${nd}.log
  echo "BASE_START m=$mm n=$nn $(date '+%T')" >> "$MASTER"
  timeout 5400 nice -n 10 magma -b Mnum:=$mn Mden:=$md Nnum:=$nnu Nden:=$nd \
      code/claude_prod_09_88_j1exact.m > "$log" 2>&1
  rc=$?
  if [ $rc -eq 124 ]; then
    echo "BASE_TIMEOUT m=$mm n=$nn (log likely empty: -b buffering)" >> "$MASTER"
  else
    up=$(grep -c "J1 UPGRADE" "$log" 2>/dev/null)
    echo "BASE_DONE m=$mm n=$nn rc=$rc upgrades=${up:-0}" >> "$MASTER"
  fi
done
echo "J1EXACT_SWEEP_DONE $(date '+%F %T')" >> "$MASTER"
