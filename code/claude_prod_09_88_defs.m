// (8,8) lane: shared Lambda_334 definitions.  Loaded (repo-root cwd) by
// claude_prod_09_88_family.m (the validation gate), j1exact, liftlocus, etc.
// Formulas: Nicholls Prop 5.9.6 / Example 5.9.4 via notes/claude_top10_08_88.md;
// validated against 15 recorded ground-truth rows by
// results/claude_prod_09_88_family_rebuild.log.  Do NOT edit semantics in
// place -- downstream scripts and the validation log depend on them.

Q_88 := Rationals();
Z_88 := Integers();
P88<x> := PolynomialRing(Q_88);

// Rational square-class representative (squarefree integer, signed).
function SqfPart(q)
    assert q ne 0;
    n := Numerator(q) * Denominator(q);
    fac := Factorization(n);
    r := Sign(n) * &*[Z_88| p[1] : p in fac | IsOdd(p[2])];
    assert IsSquare(q/r);
    return r;
end function;

// Clear denominators of y^2 = g without changing the curve over Q:
// (d*y)^2 = d^2*g with d the coefficient-denominator lcm.
function IntSextic(g)
    d := LCM([Denominator(co) : co in Coefficients(g)]);
    gi := d^2 * g;
    assert &and[Denominator(co) eq 1 : co in Coefficients(gi)];
    return gi;
end function;

// Lambda_334 member.  Returns: quintic h (J2 side, roots 0,1,a,b,c), twist
// sextic g1 (J1 side, squarefree-part-of-d2 times monic f1), and a, b, c.
function Lambda334(s, t, v)
    A := s^2 - t^4 + t^2;
    denu := -s^2*t*A*v^2 + t;
    error if denu eq 0, "degenerate: u undefined";
    u := (-s^2*A*v^2 - 2*A*v - 1) / denu;
    error if t^2 eq 1, "degenerate: t^2=1";
    denb := u^2*s^2 + 1 - t^2;
    error if denb eq 0, "degenerate: b undefined";
    a := A/(1 - t^2);
    b := A/denb;
    c := t^2;
    error if #{Q_88| 0, 1, a, b, c} ne 5, "degenerate: Weierstrass collision";
    d2 := A * (s^2*u^2 + t^4 - 2*t^2 + 1)
            * (s^4*u^2 - s^2*t^2*u^2 + s^2*u^2 - t^6 + 3*t^4 - 3*t^2 + 1);
    error if d2 eq 0, "degenerate: d2=0";
    h := x*(x-1)*(x-a)*(x-b)*(x-c);
    l1 := (-a+b+c-1)*x^2 + (2*a - 2*b*c)*x + (a*b*c - a*b - a*c + b*c);
    l2 := -x^2 + b*c;
    l3 := x^2 - a;
    ff := l1*l2*l3;
    error if Degree(ff) ne 6, "degenerate: f1 not sextic";
    f1 := ff / LeadingCoefficient(ff);
    g1 := SqfPart(d2) * f1;
    error if Discriminant(g1) eq 0 or Discriminant(h) eq 0, "degenerate: disc 0";
    return h, g1, a, b, c;
end function;

// Stage-1 base parametrization: t=(m^2+1)/(2m), s=(n^2+alpha^4)/(2n).
function StageOneST(m, n)
    t := (m^2 + 1)/(2*m);
    alpha := (m^2 - 1)/(2*m);
    assert t^2 - 1 eq alpha^2;
    s := (n^2 + alpha^4)/(2*n);
    beta := (n^2 - alpha^4)/(2*n);
    assert s^2 - alpha^4 eq beta^2;
    return s, t;
end function;

// Exact torsion invariants of Jac(y^2 = g), minimized when possible.
function ExactTorsion(g)
    C := HyperellipticCurve(IntSextic(g));
    try
        Cm := ReducedMinimalWeierstrassModel(C);
        C := SimplifiedModel(Cm);
    catch e
        ;   // fall back to the raw integral model
    end try;
    T := TorsionSubgroup(Jacobian(C));
    return Invariants(T);
end function;

// Does some order-4 torsion element double to the 2-torsion class Tcls?
function HasOrder4Over(J, Tcls)
    Tgrp, mp := TorsionSubgroup(J);
    return exists{ g : g in Tgrp | Order(g) eq 4 and 2*mp(g) eq Tcls },
           Invariants(Tgrp);
end function;
