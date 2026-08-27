#!/bin/bash
# smallmem.sh QUEUEFILE — sequentially probe fresh live T5 members rho' = q^2/(q^2-p^2)
# with sibenum (covers 3 and 4, modest depths).  QUEUEFILE lines: "p q".
SB=/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/sib_B
cd $SB
while read -r p q; do
    [ -z "$p" ] && continue
    rn=$((q*q)); rd=$((q*q - p*p))
    if [ $rd -lt 0 ]; then rn=$((-rn)); rd=$((-rd)); fi
    tag="p${p}q${q}"
    for cv in 3 4; do
        log="sm_${tag}_cv${cv}.log"
        [ -f "$log" ] && grep -q "^DONE" "$log" && continue
        nice -n 10 magma -b rn:=$rn rd:=$rd cv:=$cv nmax:=200 nm2:=30 tmax:=240 sibenum.m > "$log" 2>&1
        # surface any nondegenerate full candidate immediately
        grep "FULL CANDIDATE" "$log" | grep -v "u = 1 \*" >> $SB/SMALLMEM_CANDIDATES.txt
    done
    echo "$tag done $(date +%H:%M:%S)" >> $SB/smallmem_progress.txt
done < "$1"
echo "QUEUE $1 COMPLETE" >> $SB/smallmem_progress.txt
