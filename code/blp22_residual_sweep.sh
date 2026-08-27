#!/bin/bash
# [2,22] residual-locus sweep: solve every (r,s) base of height <= 3 over Q.
# Each base is a COMPLETE enumeration of its rational residual points
# (Res_q -> rational d-roots -> q-gcd roots).  r <-> s is a mirror symmetry:
# enumerate unordered pairs.  Usage: blp22_residual_sweep.sh SHARD NSHARDS
# One Magma at a time per shard invocation (run <= 3 shards).
#
# Run: nohup ./code/blp22_residual_sweep.sh 0 3 > /dev/null 2>&1 &  (etc.)

cd torsion_jac || exit 1
SHARD=${1:-0}; NSH=${2:-1}
MASTER=results/blp22_residual_sweep_master.log

VALS="0/1 1/1 -1/1 2/1 -2/1 3/1 -3/1 1/2 -1/2 3/2 -3/2 1/3 -1/3 2/3 -2/3 4/1 -4/1 1/4 -1/4 3/4 -3/4 4/3 -4/3 5/1 -5/1 1/5 -1/5 2/5 -2/5 3/5 -3/5 4/5 -4/5 5/2 -5/2 5/3 -5/3 5/4 -5/4"

echo "SWEEP_START shard=$SHARD/$NSH $(date '+%F %T')" >> "$MASTER"
i=-1
for rr in $VALS; do
  for ss in $VALS; do
    # unordered pairs, r != s: take rr lexicographically before ss
    [ "$rr" = "$ss" ] && continue
    [ "$(printf '%s\n%s\n' "$rr" "$ss" | sort | head -1)" != "$rr" ] && continue
    i=$((i+1))
    [ $((i % NSH)) -ne "$SHARD" ] && continue
    rn=${rr%%/*}; rd=${rr##*/}; sn=${ss%%/*}; sd=${ss##*/}
    log=results/blp22_residual_r${rn}_${rd}_s${sn}_${sd}.log
    [ -s "$log" ] && grep -q "RESIDUAL_SOLVE_DONE" "$log" && continue
    timeout 1800 nice -n 10 magma -b Rnum:=$rn Rden:=$rd Snum:=$sn Sden:=$sd \
        code/blp22_residual_solve.m > "$log" 2>&1
    rc=$?
    if [ $rc -eq 124 ]; then
      echo "BASE_TIMEOUT r=$rr s=$ss" >> "$MASTER"
    else
      hits=$(grep -c "HIT22" "$log" 2>/dev/null)
      echo "BASE_DONE r=$rr s=$ss rc=$rc hits=${hits:-0}" >> "$MASTER"
      [ "${hits:-0}" -gt 0 ] && grep "HIT22" "$log" >> "$MASTER"
    fi
  done
done
echo "SWEEP_SHARD_DONE shard=$SHARD $(date '+%F %T')" >> "$MASTER"
