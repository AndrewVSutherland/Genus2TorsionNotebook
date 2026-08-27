// claude_222_flynn43_param.m — parametrize the genus-0 Flynn (4,3) QF component,
// pull back to a 1-parameter family of [2,4]-type order-11 curves F_{t(tau)} =
// (x^2+u x+v)*Q4(x), and compute the geometry of the two [2,22]-upgrades:
//   (i) extra-root cover {Q4(r) = 0}  ([1,2,3]-type: still 2-rank 1 — for info),
//  (ii) full-split cover {Q4 = (x^2+Ax+B)(x^2+Cx+D)}  ([2,2,2]-type: 2-rank 2 = [2,22]!),
// (iii) Galois group of Q4 over Q(tau)  (D4 or smaller => rational Richelot kernels
//       along the whole family — second-generation isogeny move).
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();

A2<u,v> := AffineSpace(Q, 2);
g43 := u^4+(2*v-2)*u^3+(v^2-6*v+3)*u^2+(-8*v^2+10*v-2)*u+(-4*v^3+8*v^2-4*v+1);
C43 := Curve(A2, g43);
Cp := ProjectiveClosure(C43);
pt := Cp![-11/100, 2401/1600, 1];
Con, mp := Conic(Cp);  // genus-0 curve with a point: anticanonical; fallback below
print "conic ok";
// Prefer direct parametrization
P1 := Curve(ProjectiveSpace(Q,1));
try
  phi := Parametrization(Cp, pt, P1);
  print "parametrization found";
catch e
  print "Parametrization failed:", e`Object;
  quit;
end try;
FF<tau> := FunctionField(Q);
seq := DefiningPolynomials(phi);
uu := Evaluate(seq[1], [tau,1]) / Evaluate(seq[3], [tau,1]);
vv := Evaluate(seq[2], [tau,1]) / Evaluate(seq[3], [tau,1]);
printf "u(tau) = %o\nv(tau) = %o\n", uu, vv;

// recover t(tau): remainder coefficients of Flynn mod x^2+ux+v are quadratics in t
Pt<T> := PolynomialRing(FF);
Px<x> := PolynomialRing(Pt);
Fl := x^6+2*x^5+(2*T+3)*x^4+2*x^3+(T^2+1)*x^2+2*T*(1-T)*x+T^2;
qq, rr := Quotrem(Fl, x^2 + (Px!uu)*x + (Px!vv));
e1 := Coefficient(rr,1); e0 := Coefficient(rr,0);   // in Q(tau)[T]
gg := GCD(e1, e0);
printf "deg_T gcd = %o\n", Degree(gg);
assert Degree(gg) ge 1;
tt := -Coefficient(gg,0)/Coefficient(gg,1);         // t(tau)
printf "t(tau) = %o\n", tt;

// the family member over Q(tau)
PxF<X> := PolynomialRing(FF);
FlT := X^6+2*X^5+(2*tt+3)*X^4+2*X^3+(tt^2+1)*X^2+2*tt*(1-tt)*X+tt^2;
quo, rem := Quotrem(FlT, X^2 + FF!uu*X + FF!vv);
assert rem eq 0;
Q4 := quo;   // quartic over Q(tau)
printf "Q4 = %o\n", Q4;
print "Q4 irreducible over Q(tau):", IsIrreducible(Q4);
try
  G := GaloisGroup(Q4);
  printf "GaloisGroup(Q4/Q(tau)) order %o : %o\n", #G, GroupName(G);
catch e
  print "GaloisGroup failed:", e`Object;
end try;

// (i) extra-root cover: numerator of Q4(r) as plane curve in (tau, r)
R2<tv,rv> := PolynomialRing(Q, 2);
num := R2!0;
cfs := Coefficients(Q4);   // in Q(tau)
den := LCM([Denominator(c) : c in cfs]);
for i in [1..#cfs] do
  num +:= Evaluate(Numerator(cfs[i]*den), tv) * rv^(i-1);
end for;
A2b<t2,r2> := AffineSpace(Q,2);
Xroot := Curve(A2b, Evaluate(num, [t2, r2]));
print "extra-root cover:";
for comp in IrreducibleComponents(Xroot) do
  Z := ReducedSubscheme(comp);
  if Dimension(Z) eq 1 then
    CC := Curve(Z);
    ab := IsAbsolutelyIrreducible(CC);
    printf "  comp deg %o absirred %o", Degree(Z), ab;
    if ab then printf " GENUS %o", Genus(CC); end if;
    printf "\n";
  end if;
end for;

// (ii) full-split locus.  Q4 factors over Q(tau) as (X - rho)*C3 with rho
// rational (see factorization below).  Any factorization of Q4 into two
// RATIONAL quadratics partitions its four roots into two Galois-stable pairs;
// the pair containing the rational root rho is stable iff its other member is
// a rational root of C3.  Hence the full-split locus COINCIDES with the
// nontautological (genus-3) component of the extra-root cover above — no
// separate scheme decomposition is needed.  (An earlier version of this
// script decomposed a mis-scaled scheme here — lc^2 vs lc, caught in Codex's
// PR #8 review; its component list was invalid and is superseded by this
// argument.)
fac4 := Factorization(Q4);
printf "Q4 factorization degrees over Q(tau): %o\n", [ Degree(pe[1]) : pe in fac4 ];
for pe in fac4 do
  if Degree(pe[1]) eq 1 then
    printf "  rational root rho(tau) = %o\n", -Coefficient(pe[1],0)/Coefficient(pe[1],1);
  elif Degree(pe[1]) eq 3 then
    printf "  cubic factor C3 (irreducible: %o)\n", IsIrreducible(pe[1]);
  end if;
end for;
print "full-split locus = {C3 has a rational root} = genus-3 extra-root component (above)";

// sanity: specialize tau, check torsion [22]
for tau0 in [2, 3, 1/2, -4] do
  t0 := Evaluate(tt, tau0);
  P0<z> := PolynomialRing(Q);
  F0 := z^6+2*z^5+(2*t0+3)*z^4+2*z^3+(t0^2+1)*z^2+2*t0*(1-t0)*z+t0^2;
  if Discriminant(F0) eq 0 then continue; end if;
  d0 := LCM([Denominator(c) : c in Coefficients(F0)]);
  F0i := P0![ c*d0^2 : c in Coefficients(F0) ];
  inv := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(F0i))));
  printf "SANITY tau=%o t=%o torsion=%o\n", tau0, t0, inv;
end for;
print "FLYNN43_PARAM_DONE";
quit;
