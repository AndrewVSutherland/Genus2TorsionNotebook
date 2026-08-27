// The honest quintic component of the 4-Weierstrass locus (t1,t2)=(2,3):
// genus, rational points, and torsion validation of resulting fibers.
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

function WPdata(t)
    q := t^2/(t^2-1); m := t/(t^2-1);
    return (q-1)*(-(2*q-1)+2*m)/q, q;
end function;

// rebuild F
K<t3,t4> := RationalFunctionField(Q, 2);
t1 := K!2; t2 := K!3;
lams := []; qs := [];
for t in [t1,t2,t3,t4] do
    lam, q := WPdata(t); Append(~lams, lam); Append(~qs, q);
end for;
M := Matrix(K, 4, 4, [ [lams[i]^2, lams[i], 1, qs[i]*lams[i]] : i in [1..4] ]);
NR<a,b> := PolynomialRing(Q, 2);
num := NR!Numerator(Determinant(M));
fac := Factorization(num);
F := [ tt[1] : tt in fac | TotalDegree(tt[1]) ge 4 ][1];
printf "F = %o\n", F;

// genus via plane curve
A2 := AffineSpace(NR);
Cu := Curve(A2, F);
printf "geometric genus: %o\n", Genus(Cu);

// rational point scan: F is cubic in b for fixed a
Pb<bb> := PolynomialRing(Q);
hts := [];
for n in [-40..40] do for d in [1..8] do
    if GCD(Abs(n), d) eq 1 then Append(~hts, n/d); end if;
end for; end for;
hts := Sort(Setseq(Seqset(hts)));
pts := [];
for av in hts do
    Fb := &+[ Evaluate(Coefficient(F, b, j), [av, 0]) * bb^j : j in [0..3] ];
    if Fb eq 0 then continue; end if;
    for r in Roots(Fb) do
        bv := r[1];
        if bv in {Q!2, Q!3, av} or av in {Q!2, Q!3} then continue; end if;
        if Abs(av) in {Q!1} or Abs(bv) in {Q!1} or av eq 0 or bv eq 0 then continue; end if;
        Append(~pts, [av, bv]);
    end for;
end for;
printf "rational points found (nondeg params): %o\n", #pts;

// validate torsion on up to 8 points
nchecked := 0;
for pt in pts do
    if nchecked ge 8 then break; end if;
    ts := [Q!2, Q!3, pt[1], pt[2]];
    ls := []; qv := [];
    for t in ts do
        lam, q := WPdata(t); Append(~ls, lam); Append(~qv, q);
    end for;
    if #Set(ls) ne 4 then continue; end if;
    // quadratic through the FIRST THREE, then assert 4th automatically
    Qq := Interpolation([ls[i] : i in [1..3]], [qv[i]*ls[i] : i in [1..3]]);
    if Evaluate(Qq, ls[4]) ne qv[4]*ls[4] then continue; end if;  // det check
    g := 4*Qq^3 + (x^2-6*x+1)*Qq^2 + 2*x*(x-1)*Qq + x^2;
    if Discriminant(g) eq 0 then continue; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    C := HyperellipticCurve(den^2*g);
    tor := Invariants(TorsionSubgroup(Jacobian(C)));
    nchecked +:= 1;
    printf "  PT (%o,%o): torsion %o\n", pt[1], pt[2], tor;
end for;
printf "ROUND7B_DONE\n";
quit;
