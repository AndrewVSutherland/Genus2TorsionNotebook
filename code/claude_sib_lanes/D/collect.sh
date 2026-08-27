#!/bin/bash
# Final tally collector for Lane D.
D=/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/sib_D
echo "=== batch C-scan passes ==="
for p in A B C C2 D; do
  f=$D/batch_pass$p.log
  [ -f $f ] || continue
  last=$(grep "PASS .* DONE" $f | tail -n 1)
  if [ -n "$last" ]; then echo "pass $p: $last"; else echo "pass $p: IN FLIGHT, last: $(tail -n 1 $f)"; fi
  grep "HIT-CANDIDATE" $f | sort -u | sed 's/^/  /'
done
echo "=== sieve chunks (d>2000) ==="
cat $D/driver.log 2>/dev/null
for lg in $D/log_d*.log; do grep -h "^DONE" $lg; done 2>/dev/null
echo "=== MW probes ==="
grep -hE "member .* BND|MW-HIT" $D/mwprobe_rest.log $D/mwprobe_143_b6.log 2>/dev/null
echo "=== memberjobs ==="
echo "rd240 all-rn N=2000: $(grep -c HIT $D/rd240_allrn_n2000.txt) HIT lines, $(grep -c NEAR $D/rd240_allrn_n2000.txt) NEAR"
echo "nmdenoms odd^2 N=2000: $(grep -c . $D/nmdenoms_odd2_n2000.txt) lines"
