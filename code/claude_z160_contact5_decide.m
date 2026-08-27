// claude_z160_contact5_decide.m — FINITE decision: does any member of the reconstructed
// Elkies [32] family carry a rational 5-torsion class of the form P-infty?
// (=> torsion contains Z/32 x Z/5 = Z/160.)
// Ansatz (complete for P-infty classes on a quintic model): ftil = h5^2 - c5*(x-rho)^5
// with h5 quadratic; x^5 forces c5 = -lead(ftil). With u a VARIABLE: 5 equations in
// 5 unknowns (u, h2, h1, h0, rho) — expected dimension 0; Variety decides membership.
// (The family's marked class D has order 32, not 5, so nothing is built in.)
SetColumns(0);
SetMemoryLimit(12*10^9);
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
Q := Rationals();

// lift the integral coefficients of ftil into Q[u]
R5<uu, h2, h1, h0, rho> := PolynomialRing(Q, 5);
cf := [];
for i in [0..5] do
  ci := Coefficient(ftil, i);
  assert Denominator(ci) eq 1;
  Append(~cf, Evaluate(Numerator(ci), uu));
end for;
Pxx<X> := PolynomialRing(R5);
Fx := &+[ Pxx | cf[i+1]*X^i : i in [0..5] ];
H5 := h2*X^2 + h1*X + h0;
c5 := -cf[6];
Res := Fx - H5^2 + c5*(X - rho)^5;
eqs := [ Coefficient(Res, i) : i in [0..4] ];
assert Coefficient(Res, 5) eq 0;
I := ideal< R5 | eqs >;
// saturate away the vertical component over degenerate members (lead(ftil)(u) = 0,
// where the (x-rho)^5 term vanishes and rho becomes a free variable — this is the
// spurious dim-1 component found by the fiber probe)
print "saturating by the leading coefficient...";
t0 := Cputime();
I := Saturation(I, cf[6]);
printf "saturation done [%o s]\n", Cputime(t0);
print "computing dimension of the saturated second-5-contact system...";
t0 := Cputime();
d := Dimension(I);
printf "dimension = %o  [%o s]\n", d, Cputime(t0);
if d lt 0 then
  print "RESULT: inconsistent — NO member has a rational P-infty 5-torsion class.";
  print "THEOREM: no Z/160 upgrade of the reconstructed Elkies [32] component via P-infty 5-classes.";
elif d eq 0 then
  print "0-dimensional: computing rational points...";
  V := Variety(I);
  printf "rational solutions: %o\n", #V;
  for v in V do printf "CANDIDATE (u,h2,h1,h0,rho) = %o\n", v; end for;
  if #V eq 0 then
    print "RESULT: no rational solutions — no Z/160 member (P-infty shape) on this component.";
  end if;
else
  printf "positive dimension %o (unexpected): computing elimination to (uu,rho)\n", d;
  EI := EliminationIdeal(I, 3);
  for b in Basis(EI) do printf "ELIM: %o\n", b; end for;
end if;
print "ALL_DONE";
quit;
