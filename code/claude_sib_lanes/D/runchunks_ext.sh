#!/bin/bash
# Lane D task 1: extend the unconditional-kill enumerator beyond d=2000.
# Chunks of 100 in d; each chunk = tor22212 (-P127 cached tables) + postfilter.
# Resumable via done markers; kill by PID is safe (chunk restarts cleanly).
P=/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/prod_02_22212
D=/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/sib_D
cd $D
export OMP_NUM_THREADS=3
for lo in $(seq 2001 100 3901); do
  hi=$((lo+99))
  if [ -f $D/done_d${lo}_${hi} ]; then continue; fi
  nice -n 10 $P/tor22212 4000 $lo $hi -T $P/t127.bin -q \
      > $D/surv_d${lo}_${hi}.txt 2> $D/log_d${lo}_${hi}.log
  rc=$?
  if [ $rc -ne 0 ]; then echo "chunk $lo $hi FAILED rc=$rc" >> $D/driver.log; exit 1; fi
  nice -n 10 $P/postfilter < $D/surv_d${lo}_${hi}.txt \
      > $D/pf_d${lo}_${hi}.txt 2> $D/pf_d${lo}_${hi}.log
  touch $D/done_d${lo}_${hi}
  echo "chunk $lo $hi done: $(wc -l < $D/surv_d${lo}_${hi}.txt) surv, $(wc -l < $D/pf_d${lo}_${hi}.txt) pass pf" >> $D/driver.log
done
touch $D/ALLDONE_ext
