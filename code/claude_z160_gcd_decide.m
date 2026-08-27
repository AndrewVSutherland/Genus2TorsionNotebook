// claude_z160_gcd_decide.m — complete finite decision for the second 5-contact
// (=> Z/160) on the Elkies-32 family, via the triangular collapse:
//   h2^2 = f4 + 5 c5 rho  (x^4 eq),  h1, h0 rational in (rho, h2) with parity,
// leaving TWO univariate polynomials in rho over Q(u): EQ_A, EQ_B.
// Horizontal solutions <=> GCD nontrivial over Q(u); else candidates confined to
// rational roots of Res_rho(EQ_A, EQ_B) in u: a finite, checkable set.
SetColumns(0);
SetMemoryLimit(8*10^9);
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
u := Ku.1;
f0 := Coefficient(ftil,0); f1 := Coefficient(ftil,1); f2 := Coefficient(ftil,2);
f3 := Coefficient(ftil,3); f4 := Coefficient(ftil,4); f5 := Coefficient(ftil,5);
c5 := -f5;

// work in Q(u)[rho, H] with H = h2^2 = f4 + 5 c5 rho substituted at the end;
// track h1 = A1/(2 h2), h0 = A0/(2 h2) with A1, A0 in Q(u)[rho, H]
R2<rho, H> := PolynomialRing(Ku, 2);
FR := FieldOfFractions(R2);
Hval := FR!f4 + 5*FR!c5*rho;      // h2^2
A1 := FR!f3 + 10*FR!c5*rho^2;                       // = 2 h2 h1
h1sq := A1^2/(4*Hval);
A0 := FR!f2 - 10*FR!c5*rho^3 - h1sq;                // = 2 h2 h0
h0sq := A0^2/(4*Hval);
// x^1: f1 - 2 h1 h0 + 5 c5 rho^4 = 0 ; 2 h1 h0 = A1 A0 / (2 H)
EQA := FR!f1 - A1*A0/(2*Hval) + 5*FR!c5*rho^4;
// x^0: f0 - h0^2 - c5 rho^5 = 0
EQB := FR!f0 - h0sq - FR!c5*rho^5;
NA := Numerator(EQA); NB := Numerator(EQB);
// substitute H -> Hval is already implicit (Hval used); NA, NB should involve rho only
PuR<r1> := PolynomialRing(Ku);
hmap := hom< R2 -> PuR | [r1, PuR!0] >;   // H eliminated via Hval already
PA := hmap(NA); PB := hmap(NB);
printf "EQ_A degree in rho: %o, EQ_B degree: %o\n", Degree(PA), Degree(PB);
G := GCD(PA, PB);
printf "GCD degree in rho: %o\n", Degree(G);
if Degree(G) gt 0 then
  printf "HORIZONTAL component present: GCD = %o\n", G;
  print "factoring...";
  for t in Factorization(G) do printf "HFACTOR deg %o mult %o: %o\n", Degree(t[1]), t[2], t[1]; end for;
else
  print "no horizontal solutions: computing Res_rho for the finite candidate set...";
  RES := Resultant(PA, PB);
  num := Numerator(RES);
  printf "resultant u-degree: %o\n", Degree(num);
  rts := Roots(num);
  printf "rational u-candidates: %o\n", [r[1] : r in rts];
  Q := Rationals(); Pq<x> := PolynomialRing(Q);
  cfs := [ Numerator(Coefficient(ftil, i)) : i in [0..5] ];
  for r in rts do
    u0 := r[1];
    fv := Pq![ Evaluate(c, u0) : c in cfs ];
    if Degree(fv) lt 5 or Discriminant(fv) eq 0 then printf "u=%o DEGENERATE\n", u0; continue; end if;
    den := LCM([Denominator(c) : c in Coefficients(fv)]);
    fint := Pq![ c*den^2 : c in Coefficients(fv) ];
    Cs := SimplifiedModel(ReducedMinimalWeierstrassModel(HyperellipticCurve(fint)));
    T0 := TorsionSubgroup(Jacobian(Cs));
    printf "u=%o NONSINGULAR: torsion %o%o\n", u0, Invariants(T0),
      (#T0 mod 5 eq 0) select "  <<< 5-PART: Z/160 CANDIDATE" else "";
  end for;
end if;
print "ALL_DONE";
quit;
