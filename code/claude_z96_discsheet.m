// claude_z96_discsheet.m — resolve the disc(q3) = 0 sheet (Q1^2 = 4Q0) of the
// 3-contact system on the Elkies-32 family: doubled-support 3-torsion candidates
// f = h3^2 + kappa (x-r0)^6.  Substitute q0 = q1^2/4 (so q3 = (x+q1/2)^2), reduce
// to two equations in (E, q1) over Q(u), eliminate, and analyze.
SetColumns(0);
SetMemoryLimit(40*10^9);
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
u := Ku.1;
f0 := Coefficient(ftil,0); f1c := Coefficient(ftil,1); f2c := Coefficient(ftil,2);
f3c := Coefficient(ftil,3); f4c := Coefficient(ftil,4); f5c := Coefficient(ftil,5);

R2<E, q1> := PolynomialRing(Ku, 2);
FR := FieldOfFractions(R2);
q0 := q1^2/4;
A2 := (FR!f5c + 3*E*q1);
e2sq := A2^2/(4*E);
A1 := FR!f4c - e2sq + 3*E*(q0 + q1^2);
e1sq := A1^2/(4*E);
A0 := FR!f3c - A2*A1/(2*E) + E*(q1^3 + 6*q1*q0);
e0sq := A0^2/(4*E);
EQ1 := A2*A0/(2*E) + e1sq - E*(3*q1^2*q0 + 3*q0^2) - FR!f2c;
EQ2 := A1*A0/(2*E) - 3*E*q1*q0^2 - FR!f1c;
EQ3 := e0sq - E*q0^3 - FR!f0;
N1 := Numerator(EQ1); N2 := Numerator(EQ2); N3 := Numerator(EQ3);
printf "degrees (E,q1): %o %o %o\n", [Degree(N1,R2.i) : i in [1,2]], [Degree(N2,R2.i) : i in [1,2]], [Degree(N3,R2.i) : i in [1,2]];
t0 := Cputime();
R12 := Resultant(N1, N2, E);
R13 := Resultant(N1, N3, E);
printf "resultants done [%o s]\n", Cputime(t0);
PuR<t> := PolynomialRing(Ku);
h12 := hom< R2 -> PuR | [PuR!0, t] >;
P12 := h12(R12); P13 := h12(R13);
printf "degrees in q1: %o %o\n", Degree(P12), Degree(P13);
G := GCD(P12, P13);
printf "GCD degree in q1: %o\n", Degree(G);
if Degree(G) gt 0 then
  for tt in Factorization(G) do
    printf "DFACTOR deg %o mult %o: %o\n", Degree(tt[1]), tt[2], tt[1];
  end for;
else
  print "no horizontal disc-sheet solutions; finite u-candidates from Res:";
  RES := Resultant(P12, P13);
  num := Numerator(RES);
  printf "resultant u-degree %o\n", Degree(num);
  rts := Roots(num);
  printf "rational u-candidates: %o\n", [r[1] : r in rts];
end if;
print "ALL_DONE";
quit;
