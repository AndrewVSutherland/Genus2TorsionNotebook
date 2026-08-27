//////////////////////////////////////////////////////////////////////
// Exact ramification of the P8 map Q(u)/Q(e) at every branch value of
// the degree-12 support cover.  Run after exact reconstruction with its
// final quit omitted.
//////////////////////////////////////////////////////////////////////

print "P8_EXACT_BASE_RAMIFICATION_START";
print "P8_MAP_DEGREE",
      Max(Degree(Numerator(eP)),Degree(Denominator(eP)));

branchqs := [ Qz.1,
              Qz.1+1/3,
              Qz.1+8/15,
              Qz.1^2+2/5*Qz.1-1/15 ];
for q in branchqs do
    pull := Evaluate(q,eP);
    fn := Factorization(Numerator(pull));
    fd := Factorization(Denominator(pull));
    print "P8_PULLBACK",q,
          "ZERO_FACTORS",[<Degree(h[1]),h[2]>:h in fn],
          "POLE_FACTORS",[<Degree(h[1]),h[2]>:h in fd];
end for;

print "P8_INFINITY_PREIMAGES",
      [<Degree(h[1]),h[2]>:h in Factorization(Denominator(eP))];
print "P8_EXACT_BASE_RAMIFICATION_DONE";
quit;
