//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_target.m     lane 9 ([6,12])   2026-07-25
//
// PRIME SELECTION for the Abel-Prym MW sieve, done on the cheap side.
//
// claude_ov_m612ap_satcheck.m established that, at every prime tested,
//     ord(W mod q)  =  ord(g0 mod q)   in  Pr(F_q) = J(D)(F_q),
// where D : y^2 = -3x^6+24x^3-75 and g0 = (x^2+2x+4, 5x+5) generates J(D)(Q).
// ord(g0 mod q) costs milliseconds on the genus-2 side, while the sieve's own
// ClassGroup of the genus-4 function field E8 costs 15-60 minutes at q ~ 130
// and grows fast.  So DO NOT pick sieve primes blindly: pick the q whose
// ord(g0 mod q) is divisible by the primes ell whose admissible residue set
// the sieve has not yet collapsed to {0}.
//
// Set WANT to those ell (read them off claude_ov_m612ap_combine2.py's output:
// the rows with "ell^0 | d") and run.
//
// Usage: code/claude_magma_slot.sh -b qlo:=107 qhi:=400 code/claude_ov_m612ap_target.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetMemoryLimit(4*10^9);
if not assigned qlo then qlo := 107; elif Type(qlo) eq MonStgElt then qlo := StringToInteger(qlo); end if;
if not assigned qhi then qhi := 400; elif Type(qhi) eq MonStgElt then qhi := StringToInteger(qhi); end if;

Q := Rationals(); Px<x> := PolynomialRing(Q);
fD := -3*x^6 + 24*x^3 - 75;
JD := Jacobian(HyperellipticCurve(fD));
g0 := [P : P in Points(JD : Bound := 2000) | Order(P) eq 0][1];
WANT := {11, 13, 17, 31, 37, 43};
printf "D : y^2 = %o,  g0 = (%o, %o)\n", fD, g0[1], g0[2];
printf "looking for q in [%o..%o] whose ord(g0 mod q) is divisible by one of %o\n",
   qlo, qhi, Sort(Setseq(WANT));
print "q  ord(g0 mod q)  factorisation  wanted-primes-hit";
for q in PrimesInInterval(qlo, qhi) do
  kq := GF(q); Pq := PolynomialRing(kq);
  fq := Pq![kq!c : c in Coefficients(fD)];
  if Discriminant(fq) eq 0 then continue; end if;
  Jq := Jacobian(HyperellipticCurve(fq));
  aq := Pq![kq!c : c in Coefficients(g0[1])];
  bq := Pq![kq!c : c in Coefficients(g0[2])];
  o := Order(elt< Jq | aq, bq >);
  hit := { p[1] : p in Factorization(o) } meet WANT;
  if #hit gt 0 then
    printf "%-5o %-12o %-26o %o\n", q, o, Factorization(o), Sort(Setseq(hit));
  end if;
end for;
print "TARGET_DONE";
quit;
