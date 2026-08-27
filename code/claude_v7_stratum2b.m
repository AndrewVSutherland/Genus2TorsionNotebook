// claude_v7_stratum2b.m — stage B on the E0-hypersurface: continue the CF of
// sqrt(f) over the function field of X0 = V(E0) (where deg Q_2 = 1 holds
// identically), with pattern [3,1,2,1]; the closure condition (Q_4 constant)
// cuts out the anchor stratum S of the order-7 locus.  Output: S's second
// equation as a polynomial (numerator) pulled back to Q[a,b,c,d], written to
// data/ for the geometry pass.
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();
R4<aa,bb,cc,dd4> := PolynomialRing(Q, 4);
load "data/claude_v7_stratum2_eqs.txt";  // provides E0
// function field of V(E0) as Q(a,b,c)[d]/(E0): E0 is irreducible, deg_d = 5
F3<a,b,c> := FunctionField(Q, 3);
Pd<dvar> := PolynomialRing(F3);
E0d := &+[ Evaluate(Coefficient(E0, 4, k), [a,b,c,0]) * dvar^k : k in [0..5] ];
assert IsIrreducible(E0d);
FX<d> := quo< Pd | E0d >;
P<x> := PolynomialRing(FX);

f := x*(x-1)*(x-FX!a)*(x-FX!b)*(x^2+(FX!c)*x+d);
s0 := x^3;
for k in [1..3] do
  dd := f - s0^2;
  if Degree(dd) le 2 then break; end if;
  s0 := s0 + (Coefficient(dd, 6-k)/(2*Coefficient(s0,3)))*x^(3-k);
end for;

Pi := P!0; Qi := P!1;
pat := [];
for i in [0..3] do
  // on X0 the degree of Q_2 drops to 1; Magma's div handles it as long as
  // leading coefficients normalize in FX
  ai := (Pi + s0) div Qi;
  Append(~pat, Degree(ai));
  printf "step %o: deg a=%o deg Q=%o\n", i, Degree(ai), Degree(Qi);
  Pn := ai*Qi - Pi;
  rem := (f - Pn^2) mod Qi;
  assert rem eq 0;
  Qn := (f - Pn^2) div Qi;
  Pi := Pn; Qi := Qn;
end for;
printf "pattern %o ; deg Q_4 = %o\n", pat, Degree(Qi);
// closure: all positive-degree coefficients of Q_4 must vanish on S.
// Each coefficient is an element of FX = Q(a,b,c)[d]/(E0): lift to Pd,
// clear denominators, and map back to R4 as the closure polynomial.
h34 := hom< F3 -> FieldOfFractions(R4) | [aa,bb,cc] >;
eqs := [];
for j in [1..Degree(Qi)] do
  cj := Coefficient(Qi, j);
  if cj ne 0 then
    lft := Pd!cj;   // representative of degree <= 4 in d
    den := LCM([Denominator(Coefficient(lft,k)) : k in [0..Degree(lft)]]);
    nu := R4!0;
    for k in [0..Degree(lft)] do
      ck := Coefficient(lft, k) * den;   // in F3, now polynomial
      nu +:= R4!(Numerator(h34(ck))) * dd4^k;
    end for;
    Append(~eqs, nu);
  end if;
end for;
printf "closure equations: %o\n", #eqs;
for e in eqs do
  printf "  terms %o degrees (%o,%o,%o,%o)\n", #Terms(e), Degree(e,1), Degree(e,2), Degree(e,3), Degree(e,4);
end for;
out := Sprintf("E0 := %o;\n", E0);
for i in [1..#eqs] do out cat:= Sprintf("EC%o := %o;\n", i, eqs[i]); end for;
Write("data/claude_v7_stratum2b_eqs.txt", out : Overwrite := true);

// verify anchors on ALL closure equations
anchors := [
 [-3/5,16/15,-4/15,-52/75], [-5/3,-16/9,4/9,-52/27], [15/16,-9/16,-1/4,-39/64],
 [8/5,-1/15,-26/15,1/25], [5/8,-1/24,-13/12,1/64], [-15,-24,26,9],
 [8/3,25/9,-22/9,-13/27], [3/8,25/24,-11/12,-13/192], [9/25,24/25,-22/25,-39/625],
 [1/16,25/16,-7/4,9/64], [16,25,-28,36], [16/25,1/25,-28/25,36/625] ];
for A in anchors do
  vals := [ Evaluate(e, A) eq 0 : e in eqs ];
  printf "anchor %o: closure %o\n", A, vals;
end for;
print "STRATUM2B_DONE";
quit;
