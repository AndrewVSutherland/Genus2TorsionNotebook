#!/bin/bash
# usage: jackpot.sh P Q RN RD  — full verification + G2 dedupe of a candidate
P=/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/prod_02_22212
D=/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/sib_D
cd $P
nice -n 10 magma -b cls:=5 pp:=$1 qq:=$2 rn:=$3 rd:=$4 hitverify.m 2>&1 | tee $D/jackpot_$1_$2_$3_$4.log
cat > $D/g2q_$1_$2_$3_$4.m <<EOM
P<x> := PolynomialRing(Rationals());
A := [1,1,1,2,2];
u := $1/$2; r := $3/$4;
s := u*(r-1); m := 4*r*u*(u-1); n := r-1;
den := LCM([Denominator(z) : z in [s,m,n]]);
si := Integers()!(s*den); mi := Integers()!(m*den); ni := Integers()!(n*den);
g := GCD([si,mi,ni]); si div:= g; mi div:= g; ni div:= g;
B := [2*si^2-si*ni, 2*si^2+si*mi-2*si*ni-mi*ni, 2*si^2+si*mi-si*ni-mi*ni, -mi*ni, 4*si^2-4*si*ni-mi*ni];
f := &*[A[i] + B[i]*x : i in [1..5]];
C := HyperellipticCurve(f);
B0 := [282322361376, -8243383980, -64241207724, -114724491840, 561915878400];
f0 := &*[A[i] + B0[i]*x : i in [1..5]];
printf "SAME CURVE AS KNOWN HIT: %o\n", G2Invariants(C) eq G2Invariants(HyperellipticCurve(f0));
quit;
EOM
nice -n 10 magma -b $D/g2q_$1_$2_$3_$4.m 2>&1 | tee -a $D/jackpot_$1_$2_$3_$4.log
