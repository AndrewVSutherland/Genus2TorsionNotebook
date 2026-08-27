#!/bin/bash
# tasked_rest.sh — remaining tasked member/covers, sequential in one slot:
# j7: -25/551 cover 4 (rank 2), j8: -169/1431 cover 4 (rank 1),
# j9: 841/697 cover 4 (rank 1), j10: -169/1431 cover 3 (rank 2, smaller box),
# j11: 841/697 cover 3 (rank 3, N3 box).
SB=/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/sib_B
cd $SB
nice -n 10 magma -b rn:=-25 rd:=551 cv:=4 nmax:=300 nm2:=60 tmax:=800 sibenum.m > j7_m25_551_cv4.log 2>&1
nice -n 10 magma -b rn:=-169 rd:=1431 cv:=4 nmax:=300 nm2:=60 tmax:=500 sibenum.m > j8_m169_1431_cv4.log 2>&1
nice -n 10 magma -b rn:=841 rd:=697 cv:=4 nmax:=300 nm2:=60 tmax:=500 sibenum.m > j9_m841_697_cv4.log 2>&1
nice -n 10 magma -b rn:=-169 rd:=1431 cv:=3 nmax:=300 nm2:=50 tmax:=700 sibenum.m > j10_m169_1431_cv3.log 2>&1
nice -n 10 magma -b rn:=841 rd:=697 cv:=3 nmax:=300 nm2:=14 tmax:=700 sibenum.m > j11_m841_697_cv3.log 2>&1
echo "TASKED_REST COMPLETE" >> $SB/smallmem_progress.txt
