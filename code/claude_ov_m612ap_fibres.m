//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_fibres.m     lane 9 ([6,12])   2026-07-25
//
// The two special fibres of E8 -> P^1_x over Q, which the height-bounded
// rational-point scan (code/claude_ov_m612ap_ptscan.gp) cannot see:
//   x = 0        the boundary, where the two iota-fixed rational points live;
//   x = infinity not covered by a scan over x0 = a/b at all.
//
// Result:
//   over x = infinity : ONE place, degree 8  =>  E8 has NO rational point
//                       over x = infinity;
//   over x = 0        : TWO places, both degree 1, ramification 2 and 6
//                       =>  exactly the two boundary points b1, b2, and this
//                       confirms over Q the (2,6) pattern that the sieve
//                       asserts mod q, plus r = 2 from Riemann-Hurwitz.
//
// Usage: code/claude_magma_slot.sh -b code/claude_ov_m612ap_fibres.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetMemoryLimit(6*10^9);
QQ := Rationals();
Fx<X> := RationalFunctionField(QQ);
Py<Y> := PolynomialRing(Fx);
Q8 := Y^8 + (216*X^4+72*X^3-24*X^2)*Y^4
    + (-1296*X^6-1728*X^5-432*X^4+64*X^3)*Y^2
    + (-3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4);
F<yy> := FunctionField(Q8);
printf "genus %o\n", Genus(F);

pinf := Poles(F!X);
printf "places over x = infinity: %o\n", #pinf;
for P in pinf do
  printf "  degree %o, ramification %o\n", Degree(P), RamificationIndex(P);
end for;
printf "rational points over x = infinity: %o\n", #[P : P in pinf | Degree(P) eq 1];

z0 := Zeros(Fx.1)[1];
d0 := Decomposition(F, z0);
d0 := [ Type(D) eq Tup select D[1] else D : D in d0 ];
printf "places over x = 0: %o\n", #d0;
for P in d0 do
  printf "  degree %o, ramification %o\n", Degree(P), RamificationIndex(P);
end for;
printf "rational points over x = 0: %o\n", #[P : P in d0 | Degree(P) eq 1];
print "FIBRES_DONE";
quit;
