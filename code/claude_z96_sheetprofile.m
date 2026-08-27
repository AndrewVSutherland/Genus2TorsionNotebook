// claude_z96_sheetprofile.m — factor-degree profile of the 3-contact cover of the
// Elkies-32 family WITHOUT the huge bivariate factorization: specialize u = u0,
// redo the (cheap) elimination over Q, factor the degree-~160 univariate eliminant,
// and compare degree multisets across several u0.  Stable small degrees across
// specializations = candidate low-degree sheets (worth the full computation);
// none = the cover has only high-degree sheets (Faltings-finite, structurally closed).
SetColumns(0);
SetMemoryLimit(8*10^9);
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
Q := Rationals();
cfs := [ Numerator(Coefficient(ftil, i)) : i in [0..5] ];

for u0 in [Q| 5, -7/2, 11/3, 13, -17/5, 23/2 ] do
  fv := [ Evaluate(c, u0) : c in cfs ];
  if fv[6] eq 0 then continue; end if;
  Pq<x> := PolynomialRing(Q);
  fpol := Pq!fv;
  if Discriminant(fpol) eq 0 then continue; end if;
  R3<E, q1, q0> := PolynomialRing(Q, 3);
  FR3 := FieldOfFractions(R3);
  f5 := FR3!fv[6]; f4 := FR3!fv[5]; f3 := FR3!fv[4]; f2 := FR3!fv[3]; f1 := FR3!fv[2]; f0 := FR3!fv[1];
  A2 := (f5 + 3*E*q1);
  e2sq := A2^2/(4*E);
  A1 := f4 - e2sq + 3*E*(q0 + q1^2);
  e1sq := A1^2/(4*E);
  A0 := f3 - A2*A1/(2*E) + E*(q1^3 + 6*q1*q0);
  e0sq := A0^2/(4*E);
  EQ1 := A2*A0/(2*E) + e1sq - E*(3*q1^2*q0 + 3*q0^2) - f2;
  EQ2 := A1*A0/(2*E) - 3*E*q1*q0^2 - f1;
  EQ3 := e0sq - E*q0^3 - f0;
  N1 := Numerator(EQ1); N2 := Numerator(EQ2); N3 := Numerator(EQ3);
  R12 := Resultant(N1, N2, E);
  R13 := Resultant(N1, N3, E);
  G := GCD(R12, R13);
  A := R12 div G; B := R13 div G;
  RF := Resultant(A, B, q0);
  if RF eq 0 then printf "u0=%o: cofactor resultant zero (deeper common factor)\n", u0; continue; end if;
  Pq1<T> := PolynomialRing(Q);
  hm := hom< R3 -> Pq1 | [Pq1!0, T, Pq1!0] >;
  RFq := hm(RF);
  sq := GCD(RFq, Derivative(RFq));
  RFs := RFq div sq;
  fac := Factorization(RFs);
  printf "u0=%o: eliminant deg %o, squarefree %o, factor degrees %o\n",
    u0, Degree(RFq), Degree(RFs), Sort([Degree(t[1]) : t in fac]);
  // also profile the GCD part (shared locus) at this fiber
  Gq := hm(Resultant(G, R3.3 - 0, q0));  // G may involve q0; project crudely
  if Degree(G, q0) eq 0 and Degree(G, q1) gt 0 then
    Gq1 := hm(G);
    facG := Factorization(Gq1 div GCD(Gq1, Derivative(Gq1)));
    printf "u0=%o: G-part factor degrees %o\n", u0, Sort([Degree(t[1]) : t in facG]);
  end if;
end for;
print "ALL_DONE";
quit;
