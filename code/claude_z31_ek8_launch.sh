#!/bin/bash
# claude_z31_ek8_launch.sh — production driver for the EK disc-8 31-torsion
# sieve on aws-spot-11 (192 cores).  NOT launched by the build session.
#
# Two-pass design over the H=60 slice box (~4378 r-slices):
#   pass 1: 4 CRT primes, recbound 7e5  (unique-reconstruction regime, fast)
#   pass 2: 5 CRT primes, recbound 2e6  (deep pass, ~11 h wall)
# Both passes re-run with a DISJOINT second prime set to cover solved s with
# denominators divisible by a pass-A prime.  Shards of 500 slices run
# sequentially (one process uses all cores); logs are append-safe so a spot
# interruption resumes at the first shard whose log lacks SWEEP_DONE.
# Exact stage afterwards:
#   grep -h "^CAND" ek8_prod_*.log | sort -u  ->  per line run
#   magma -b Rn:=<rn> Rd:=<rd> Sn:=<sn> Sd:=<sd> claude_z31_ek8_exact.m
set -e
cd ~/z31
gcc -O3 -march=native -fopenmp -o ek8_sieve claude_z31_ek8_sieve.c
./ek8_sieve selftest | tail -1 | grep -q "fail=0" || { echo SELFTEST_FAILED; exit 1; }
H=60
NSL=4400                      # upper bound on slice count; sweep clips to #slices
STEP=500
# prime set A (all = 3 mod 4, > 5, != 31)
A_CRT="1019 1031 1039 1051"; A_VER="1063 1087 1103 1123 1151 1163 1187 1223"
# prime set B (disjoint)
B_CRT="1259 1283 1291 1303"; B_VER="1307 1319 1327 1367 1399 1423 1427 1439"
run_pass () { # $1=tag $2=ncrt $3=crt $4=nver $5=ver $6=recbound
  local tag=$1 ncrt=$2 crt="$3" nver=$4 ver="$5" rb=$6
  for lo in $(seq 0 $STEP $((NSL-1))); do
    hi=$((lo+STEP))
    log=ek8_prod_${tag}_${lo}_${hi}.log
    if grep -qs SWEEP_DONE $log; then echo "skip $log (done)"; continue; fi
    echo "=== $tag shard [$lo,$hi) $(date -u +%FT%TZ)"
    ./ek8_sieve sweep $H $ncrt $crt $nver $ver $lo $hi $rb > $log 2>&1
  done
}
run_pass 1A 4 "$A_CRT" 8 "$A_VER" 700000
run_pass 1B 4 "$B_CRT" 8 "$B_VER" 700000
run_pass 2A 5 "1019 1031 1039 1051 1063" 7 "1087 1103 1123 1151 1163 1187 1223" 2000000
run_pass 2B 5 "1259 1283 1291 1303 1307" 7 "1319 1327 1367 1399 1423 1427 1439" 2000000
echo ALL_PASSES_DONE
