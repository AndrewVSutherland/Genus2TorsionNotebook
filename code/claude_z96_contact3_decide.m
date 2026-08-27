// claude_z96_contact3_decide.m — decision analysis: does any member of the reconstructed
// Elkies [32] family carry a rational 3-torsion class?  (=> cyclic Z/96)
// Every 3-torsion class on a quintic model has Mumford degree exactly 2, so the ansatz
//   ftil = h3^2 + kappa*q3^3,  q3 = x^2+q1*x+q0 monic,  kappa = -e3^2  (from x^6)
// is COMPLETE.  Triangular elimination of e2,e1,e0 leaves EQ1,EQ2,EQ3 in (E=e3^2, q1, q0)
// over Q(u); we eliminate to a plane curve and analyze its components.
SetColumns(0);
SetMemoryLimit(24*10^9);
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
u := Ku.1; x := Px.1;
printf "family built; ftil coefficient u-degrees: %o\n",
  [<i, Degree(Numerator(Coefficient(ftil,i)))> : i in [0..5]];

// quick member sanity: order-32 class present at two random parameters
Q := Rationals(); Pq<t> := PolynomialRing(Q);
for uv in [Q| 2, -3/2] do
  fv := Pq![Evaluate(Coefficient(ftil, i), uv) : i in [0..5]];
  if Degree(fv) ne 5 or Discriminant(fv) eq 0 then printf "u=%o degenerate member\n", uv; continue; end if;
  den := LCM([Denominator(c) : c in Coefficients(fv)]);
  fint := Pq![ c*den^2 : c in Coefficients(fv) ];
  T := TorsionSubgroup(Jacobian(HyperellipticCurve(fint)));
  printf "member u=%o: torsion %o\n", uv, Invariants(T);
end for;

f5 := Coefficient(ftil, 5); f4 := Coefficient(ftil, 4); f3 := Coefficient(ftil, 3);
f2 := Coefficient(ftil, 2); f1 := Coefficient(ftil, 1); f0 := Coefficient(ftil, 0);

// work over fraction field in (E, q1, q0)
R3<E, q1, q0> := PolynomialRing(Ku, 3);
FR3 := FieldOfFractions(R3);
e3sq := FR3!E;
// triangular solve (e2,e1,e0 scaled by e3: even-power bookkeeping via E)
// e2 = (f5 + 3E q1)/(2 e3); e2^2 = (f5+3Eq1)^2/(4E)
A2 := (FR3!f5 + 3*E*q1);           // = 2 e3 e2
e2sq := A2^2/(4*e3sq);
// x^4: 2 e3 e1 + e2^2 + kappa(3q0+3q1^2) = f4 ; kappa = -E
A1 := FR3!f4 - e2sq + 3*E*(q0 + q1^2);     // = 2 e3 e1
e1sq := A1^2/(4*e3sq);
// x^3: 2 e3 e0 + 2 e2 e1 + kappa(q1^3+6q1q0) = f3 ; 2 e2 e1 = 2*(A2/2e3)(A1/2e3) = A2 A1/(2E)
A0 := FR3!f3 - A2*A1/(2*e3sq) + E*(q1^3 + 6*q1*q0);   // = 2 e3 e0
e0sq := A0^2/(4*e3sq);
// EQ1 (x^2): 2 e2 e0 + e1^2 - E(3q1^2 q0 + 3 q0^2) - f2 = 0 ; 2 e2 e0 = A2 A0/(2E)
EQ1 := A2*A0/(2*e3sq) + e1sq - E*(3*q1^2*q0 + 3*q0^2) - FR3!f2;
// EQ2 (x^1): 2 e1 e0 - 3E q1 q0^2 - f1 = 0 ; 2 e1 e0 = A1 A0/(2E)
EQ2 := A1*A0/(2*e3sq) - 3*E*q1*q0^2 - FR3!f1;
// EQ3 (x^0): e0^2 - E q0^3 - f0 = 0
EQ3 := e0sq - E*q0^3 - FR3!f0;

N1 := Numerator(EQ1); N2 := Numerator(EQ2); N3 := Numerator(EQ3);
printf "numerator degrees (E,q1,q0): %o %o %o\n",
  [Degree(N1, R3.i) : i in [1..3]], [Degree(N2, R3.i) : i in [1..3]], [Degree(N3, R3.i) : i in [1..3]];

// lift to Q[u, E, q1, q0] and work over Q (robust GCD/factorization)
P4<U, EE, Q1, Q0> := PolynomialRing(Rationals(), 4);
function LiftP4(N)
  res := P4!0;
  mons := Monomials(N); cfs := Coefficients(N);
  den := LCM([Denominator(c) : c in cfs]);
  for i in [1..#mons] do
    c := cfs[i]*den;
    assert Denominator(c) eq 1;
    m := mons[i];
    res +:= Evaluate(Numerator(c), U) * EE^Degree(m, R3.1) * Q1^Degree(m, R3.2) * Q0^Degree(m, R3.3);
  end for;
  return res;
end function;
L1 := LiftP4(N1); L2 := LiftP4(N2); L3 := LiftP4(N3);

t0 := Cputime();
R12 := Resultant(L1, L2, EE);
printf "Res_E(L1,L2) done [%o s], degs (U,Q1,Q0) = %o\n", Cputime(t0), [Degree(R12, P4.i) : i in [1,3,4]];
t0 := Cputime();
R13 := Resultant(L1, L3, EE);
printf "Res_E(L1,L3) done [%o s], degs (U,Q1,Q0) = %o\n", Cputime(t0), [Degree(R13, P4.i) : i in [1,3,4]];

t0 := Cputime();
G := GCD(R12, R13);
printf "GCD done [%o s], degs (U,Q1,Q0) = %o\n", Cputime(t0), [Degree(G, P4.i) : i in [1,3,4]];
if Degree(G, Q1) gt 0 or Degree(G, Q0) gt 0 then
  print "COMMON FACTOR (the shared locus) — factoring:";
  fac := Factorization(G);
  for i in [1..#fac] do
    printf "GFACTOR %o degs(U,Q1,Q0)=%o mult=%o : %o\n", i,
      [Degree(fac[i][1], P4.j) : j in [1,3,4]], fac[i][2],
      #Monomials(fac[i][1]) le 40 select fac[i][1] else "large";
  end for;
end if;
A := R12 div G; B := R13 div G;
printf "cofactor degs: A %o, B %o\n", [Degree(A, P4.i) : i in [1,3,4]], [Degree(B, P4.i) : i in [1,3,4]];
if Degree(A, Q0) gt 0 and Degree(B, Q0) gt 0 then
  t0 := Cputime();
  RF := Resultant(A, B, Q0);
  printf "Res_q0(cofactors) done [%o s], degs (U,Q1) = %o\n", Cputime(t0), [Degree(RF, P4.i) : i in [1,3]];
  if RF ne 0 then
    fac2 := Factorization(RF);
    printf "cofactor eliminant: %o factors\n", #fac2;
    for i in [1..#fac2] do
      printf "CFACTOR %o degs(U,Q1)=%o mult=%o : %o\n", i,
        [Degree(fac2[i][1], P4.j) : j in [1,3]], fac2[i][2],
        #Monomials(fac2[i][1]) le 40 select fac2[i][1] else "large";
    end for;
  else
    print "cofactor resultant ALSO zero (nested common factor)";
  end if;
end if;
print "ALL_DONE";
quit;
