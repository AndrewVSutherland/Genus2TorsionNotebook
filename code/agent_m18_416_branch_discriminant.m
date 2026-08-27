
//////////////////////////////////////////////////////////////////////
//  Branch discriminant of the [4,16] surface Sigma' : E1 = E0 = 0 in
//  A^4_{R,w,a,b} at the <0,0> double branch.
//
//  agent_m18_416_e1e0_reduction.m showed both cores have lowest form
//  b^2 (times units/R) at the origin: the surface is locally a DOUBLE
//  branch in b.  This script:
//   1. expands the cores in b: E = e0 + e1*b + e2*b^2 + e3*b^3 (+e4*b^4)
//      and reports the vanishing orders of the e_j at the origin;
//   2. forms the naive branch discriminant Delta = e1^2 - 4*e0*e2 of the
//      b-quadratic part for both cores and factors them: a square factor
//      structure here is the explicit "square root H" behind the 7-adic
//      never-smooth tower;
//   3. computes Res_b(E1core, E0core) modulo a few primes to learn the
//      factor-degree profile of the projected surface in A^3_{R,w,a}
//      (the full rational resultant may be infeasible; mod-p first).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();
R4<R,w,a,b> := PolynomialRing(Q, 4, "grevlex");
K4 := FieldOfFractions(R4);
PX<x> := PolynomialRing(K4);

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
c4 := R + 2 + 4*t;
q := x^2 + a*x + b;
D := f - c4*(x+R)*q^2;
d4 := Coefficient(D,4); d3 := Coefficient(D,3); d2 := Coefficient(D,2);
d1 := Coefficient(D,1); d0 := Coefficient(D,0);
E1 := 8*d4^2*d1 - d3*(4*d4*d2 - d3^2);
E0 := 64*d4^3*d0 - (4*d4*d2 - d3^2)^2;

function ClearPrim(g)
    gn := Numerator(g);
    den := LCM([Denominator(co) : co in Coefficients(gn)]);
    gn := den*gn;
    gg := GCD([Z!co : co in Coefficients(gn)]);
    return R4!(gn/gg);
end function;

E1n := ClearPrim(E1); E0n := ClearPrim(E0);
// strip the R-1 boundary factors found before
while IsDivisibleBy(E1n, R-1) do E1n := ExactQuotient(E1n, R-1); end while;
while IsDivisibleBy(E0n, R-1) do E0n := ExactQuotient(E0n, R-1); end while;
printf "E1core deg %o terms %o | E0core deg %o terms %o\n",
    TotalDegree(E1n), #Terms(E1n), TotalDegree(E0n), #Terms(E0n);

// ---- 1. b-expansion ----
R3<R_, w_, a_> := PolynomialRing(Q, 3, "grevlex");
h34 := hom<R4 -> PolynomialRing(R3) |
    [PolynomialRing(R3)!R_, PolynomialRing(R3)!w_, PolynomialRing(R3)!a_,
     PolynomialRing(R3).1]>;
E1b := h34(E1n);   // polynomial in b over Q[R,w,a]
E0b := h34(E0n);
printf "deg_b E1core = %o, deg_b E0core = %o\n", Degree(E1b), Degree(E0b);

function OrdAtOrigin(g)   // g in R3
    if g eq 0 then return -1; end if;
    return Min([TotalDegree(m) : m in Terms(g)]);
end function;

print "b-coefficient vanishing orders at origin (E1core):";
for j in [0..Degree(E1b)] do
    cj := Coefficient(E1b, j);
    printf "  e1_%o: deg %o, ord_origin %o, terms %o\n",
        j, TotalDegree(cj), OrdAtOrigin(cj), #Terms(cj);
end for;
print "b-coefficient vanishing orders at origin (E0core):";
for j in [0..Degree(E0b)] do
    cj := Coefficient(E0b, j);
    printf "  e0_%o: deg %o, ord_origin %o, terms %o\n",
        j, TotalDegree(cj), OrdAtOrigin(cj), #Terms(cj);
end for;

// ---- 2. branch discriminants of the b-quadratic parts ----
for pair in [<"E1core", E1b>, <"E0core", E0b>] do
    name := pair[1]; Eb := pair[2];
    ee0 := Coefficient(Eb, 0); ee1 := Coefficient(Eb, 1); ee2 := Coefficient(Eb, 2);
    Delta := ee1^2 - 4*ee0*ee2;
    if Delta eq 0 then
        printf "%o: Delta == 0 (perfect square quadratic!)\n", name;
        continue;
    end if;
    printf "%o: Delta deg %o terms %o\n", name, TotalDegree(Delta), #Terms(Delta);
    facD := Factorization(Delta);
    printf "%o Delta factorization degrees/mults:\n", name;
    for ffd in facD do
        printf "  <deg %o, terms %o, mult %o>", TotalDegree(ffd[1]), #Terms(ffd[1]), ffd[2];
        if #Terms(ffd[1]) le 10 then printf "  %o", ffd[1]; end if;
        printf "\n";
    end for;
    sqfree := &*([R3!1] cat [ffd[1] : ffd in facD | IsOdd(ffd[2])]);
    printf "%o Delta squarefree part: deg %o terms %o%o\n",
        name, TotalDegree(sqfree), #Terms(sqfree),
        (#Terms(sqfree) le 12 select Sprintf("  %o", sqfree) else "");
end for;

// ---- 3. mod-p resultant profile of Res_b(E1core, E0core) ----
for p in [101, 103] do
    Fp := GF(p);
    R3p<Rp, wp, ap> := PolynomialRing(Fp, 3, "grevlex");
    Pb := PolynomialRing(R3p);
    hp := hom<R4 -> Pb | [Pb!Rp, Pb!wp, Pb!ap, Pb.1]>;
    E1p := hp(E1n); E0p := hp(E0n);
    tim := Cputime();
    Res := Resultant(E1p, E0p);
    printf "p=%o: Res_b degree %o, terms %o (%.1o s)\n",
        p, TotalDegree(Res), #Terms(Res), Cputime(tim);
    tim := Cputime();
    facR := Factorization(Res);
    printf "p=%o: Res_b factor profile (%.1o s):\n", p, Cputime(tim);
    for fr in facR do
        printf "  <deg %o, mult %o>", TotalDegree(fr[1]), fr[2];
        if #Terms(fr[1]) le 8 then printf "  %o", fr[1]; end if;
        printf "\n";
    end for;
end for;

print "DONE";
quit;
