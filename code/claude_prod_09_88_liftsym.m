// (8,8) lane, step 1.1 part B: SYMBOLIC order-4 sections over the class
// T_j on the full Lambda_334 chart Q(s,t,v), via the algebraic halving
// identity (cf. the halving-and-doubling machinery):
//   D in J1[4] with 2D = T_{lj}   <==>   exists cubic beta, quadratic u, kappa:
//       beta^2 - g1 = kappa * u^2 * lj ,   g1 = gamma * l1*l2*l3
// Evaluating at a root of lj forces lj | beta, i.e. beta = lj*(mu*x+nu), and
// dividing by lj the identity becomes: the quartic
//       Q4(x; mu,nu) := lj*(mu*x+nu)^2 - gamma * (product of the other two l's)
// must be kappa * (monic quadratic)^2.  Matching coefficients gives kappa=e4,
// p = e3/(2 e4), q = (4 e4 e2 - e3^2)/(8 e4^2) and TWO closure equations
//   F1 := 8 e4^2 e1 - e3 (4 e4 e2 - e3^2) = 0
//   F2 := 64 e4^3 e0 - (4 e4 e2 - e3^2)^2 = 0
// in (mu,nu) over K = Q(s,t,v).  Rational solutions <=> rational sections;
// u_D = x^2 + p x + q is the Mumford u-polynomial of the section.
// The class SOLVED here is T3 = [l3] = [x^2 - a]: the numeric survey
// (results/claude_prod_09_88_liftlocus_partA.log) shows the norm-presieve
// passers at stage-1 members all live over l3.
//
// VALIDATION: the symbolic sections are specialized at member
// (m,n,v) = (2,1,1) and compared against the exact order-4 torsion of J1.
//
// Output: the section u-polynomials over K (printed), for the lambda-cover
// construction and genus gate (step 1.2).
//
// Run: nohup magma -b code/claude_prod_09_88_liftsym.m > results/claude_prod_09_88_liftsym.log 2>&1 &

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

K<s,t,v> := RationalFunctionField(Rationals(), 3);

// Lambda_334 family data over K (same formulas as code/claude_prod_09_88_defs.m)
A := s^2 - t^4 + t^2;
u334 := (-s^2*A*v^2 - 2*A*v - 1) / (-s^2*t*A*v^2 + t);
a := A/(1 - t^2);
b := A/(u334^2*s^2 + 1 - t^2);
c := t^2;
d2 := A * (s^2*u334^2 + t^4 - 2*t^2 + 1)
        * (s^4*u334^2 - s^2*t^2*u334^2 + s^2*u334^2 - t^6 + 3*t^4 - 3*t^2 + 1);

PR<mu,nu> := PolynomialRing(K, 2);
Px<x> := PolynomialRing(PR);

l1 := (-a+b+c-1)*x^2 + (2*a - 2*b*c)*x + (a*b*c - a*b - a*c + b*c);
l2 := -x^2 + b*c;
l3 := x^2 - a;
lead12 := K!(LeadingCoefficient(l1)*LeadingCoefficient(l2));
// g1 = d2 * f1 = (d2/lead(l1 l2 l3)) * l1*l2*l3 ; lead(l3)=1.
gamma := d2 / lead12;

printf "CHART a=%o\nCHART c=%o\n", a, c;

// Halving of T3 = [l3]: Q4 = l3*(mu*x+nu)^2 - gamma*l1*l2
Q4 := l3*(mu*x + nu)^2 - (Px!gamma)*l1*l2;
assert Degree(Q4) le 4;
e4 := Coefficient(Q4, 4); e3 := Coefficient(Q4, 3);
e2 := Coefficient(Q4, 2); e1 := Coefficient(Q4, 1); e0 := Coefficient(Q4, 0);

F1 := 8*e4^2*e1 - e3*(4*e4*e2 - e3^2);
F2 := 64*e4^3*e0 - (4*e4*e2 - e3^2)^2;
printf "CLOSURE_DEGREES F1=%o F2=%o\n", TotalDegree(F1), TotalDegree(F2);

// Eliminate nu.
R := Resultant(F1, F2, nu);
printf "RESULTANT_COMPUTED degree_in_mu=%o\n", Degree(R, mu);

// Rational sections: rational roots of R as a polynomial in mu over K.
Pmu := PolynomialRing(K);
Rmu := &+[Pmu| K!Coefficient(R, mu, i) * Pmu.1^i : i in [0..Degree(R, mu)]];
Rmu := Rmu / LeadingCoefficient(Rmu);
fac := Factorization(Rmu);
printf "RESULTANT_FACTOR_DEGREES %o\n", [<Degree(fp[1]), fp[2]> : fp in fac];

// collect candidate rational mu values (roots of linear factors)
mucands := [ -Coefficient(fp[1], 0) : fp in fac | Degree(fp[1]) eq 1 ];
mucands := [ mv : mv in mucands | mv ne 0 ];
printf "RATIONAL_MU_CANDIDATES %o\n", #mucands;

// For each candidate mu, solve F1(mu, nu) = 0 for rational nu, verify F2,
// and reconstruct the section u_D = x^2 + p x + q.
sections := [];   // tuples <muv, nuv, p, q>
for muv in mucands do
    F1n := UnivariatePolynomial(Evaluate(F1, mu, muv));
    F2n := UnivariatePolynomial(Evaluate(F2, mu, muv));
    if F1n eq 0 then continue; end if;
    rts := [ r[1] : r in Roots(F1n) ];
    for nuv in rts do
        if Evaluate(F2n, nuv) ne 0 then continue; end if;
        // reconstruct kappa, p, q at (muv, nuv)
        ev := func< e | Evaluate(Evaluate(e, mu, muv), nu, nuv) >;
        k4 := K!ev(e4); k3 := K!ev(e3); k2 := K!ev(e2);
        if k4 eq 0 then continue; end if;
        pp := k3/(2*k4);
        qq := (4*k4*k2 - k3^2)/(8*k4^2);
        Append(~sections, <muv, nuv, pp, qq>);
        printf "SECTION mu=%o\nSECTION nu=%o\nSECTION p=%o\nSECTION q=%o\n",
               muv, nuv, pp, qq;
    end for;
end for;
printf "N_SECTIONS %o\n", #sections;
error if #sections eq 0, "no rational sections found symbolically";

// VALIDATION at member (m,n,v) = (2,1,1):
// stage-1 base: t0=(m^2+1)/(2m), alpha=(m^2-1)/(2m), s0=(n^2+alpha^4)/(2n)
m0 := 2/1; n0 := 1/1; v0 := 1/1;
t0 := (m0^2+1)/(2*m0); al := (m0^2-1)/(2*m0); s0 := (n0^2+al^4)/(2*n0);
QQ := Rationals(); PQ<xr> := PolynomialRing(QQ);
spec := func< f | Evaluate(f, [s0, t0, v0]) >;
// specialize section u-polys
uspec := { PQ | xr^2 + spec(sec[3])*xr + spec(sec[4]) : sec in sections };
printf "VALIDATION_SPECIALIZED_U %o\n", uspec;

// exact order-4 torsion over T3 at the member, via the defs machinery
load "code/claude_prod_09_88_defs.m";
sv, tv := StageOneST(m0, n0);
assert sv eq s0 and tv eq t0;
h, g1m, am, bm, cm := Lambda334(s0, t0, v0);
l3m := P88.1^2 - am;
J1 := Jacobian(HyperellipticCurve(IntSextic(g1m)));
Tg, mp := TorsionSubgroup(J1);
assert Invariants(Tg) eq [4,4];
T3pt := J1![l3m, P88!0];
utrue := { };
for g in Tg do
    if Order(g) eq 4 and 2*mp(g) eq T3pt then
        Include(~utrue, mp(g)[1]);
    end if;
end for;
printf "VALIDATION_TRUE_U %o\n", utrue;
// compare (coerce both into PQ)
uspec2 := { PQ | pol : pol in uspec };
utrue2 := { PQ![Coefficient(pol,i) : i in [0..Degree(pol)]] : pol in utrue };
match := uspec2 eq utrue2;
printf "VALIDATION_MATCH %o\n", match;
assert match;
printf "LIFTSYM_T3_SECTIONS_VALIDATED\n";
quit;
