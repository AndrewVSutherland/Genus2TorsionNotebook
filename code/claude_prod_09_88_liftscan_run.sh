#!/bin/bash
# (8,8) lane: base-grid lambda-compatibility scan (liftgate.m Scan mode).
# The lift cover has uniform v-fiber genus 3 (bases (2,1),(3,1/3),(2,3)), so
# rational lift points are sporadic per base; this sweeps many bases with the
# exact necessary condition.  Any LIFTHIT is an (8,4)-candidate -> follow up
# with exact TorsionSubgroup(J1).
# One Magma at a time (ONE job slot).  Per-base timeout 900 s.
#
# Run: nohup ./code/claude_prod_09_88_liftscan_run.sh > /dev/null 2>&1 &

cd torsion_jac || exit 1
HV=${1:-60}
MASTER=results/claude_prod_09_88_liftscan_master.log

MS="2/1 3/1 4/1 5/1 1/2 1/3 2/3 3/2 5/2 2/5 4/3 3/4"
NS="1/1 2/1 3/1 5/1 1/2 1/3 2/5 5/2"

echo "LIFTSCAN_SWEEP_START $(date '+%F %T') Hv=$HV" >> "$MASTER"
for mm in $MS; do
  for nn in $NS; do
    mn=${mm%%/*}; md=${mm##*/}
    nnu=${nn%%/*}; nd=${nn##*/}
    log=results/claude_prod_09_88_liftscan_h${HV}_m${mn}_${md}_n${nnu}_${nd}.log
    timeout $((300 + HV*10)) nice -n 10 magma -b Mnum:=$mn Mden:=$md Nnum:=$nnu Nden:=$nd \
        Lj:=3 Scan:=1 Hv:=$HV code/claude_prod_09_88_liftgate.m > "$log" 2>&1
    rc=$?
    if [ $rc -eq 124 ]; then
      echo "BASE_TIMEOUT m=$mm n=$nn" >> "$MASTER"
    else
      hits=$(grep -c "LIFTHIT" "$log" 2>/dev/null)
      echo "BASE_DONE m=$mm n=$nn rc=$rc hits=${hits:-0}" >> "$MASTER"
      if [ "${hits:-0}" -gt 0 ]; then
        grep "LIFTHIT" "$log" >> "$MASTER"
      fi
    fi
  done
done
echo "LIFTSCAN_SWEEP_DONE $(date '+%F %T')" >> "$MASTER"
