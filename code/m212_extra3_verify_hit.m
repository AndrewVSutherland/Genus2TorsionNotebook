//////////////////////////////////////////////////////////////////////
//  Verify and print an independent 3-torsion point on M(2,12).
//
//  The first hit from m212_extra3_search.m is
//      z = -5/3, r = -3/5, a = -10/9.
//
//  On the odd model obtained by moving the chosen root of T+1 to
//  infinity, the marked point gives D of order 12.  The built-in
//  3-torsion is 4D.  This script prints an order-3 class not equal to
//  0, 4D, or -4D.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<X> := PolynomialRing(Q);

function OddQuinticAtRoot(W, w)
    out := P!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            out +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
        end for;
    end for;
    return out;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

z := Q!-5/3;
r := Q!-3/5;
a := (1-z^2)/(4*(r+1));

T := a*X^2 - X + r;
h := (X-r)*(T+1);
W := h^2 + 4*a*X^2*T*(T+1);
w := 2*(r+1)/(1+z);

assert a eq Q!-10/9;
assert Evaluate(T+1, w) eq 0;

f5 := OddQuinticAtRoot(W, w);
Xp := -1/w;
Yp := Evaluate(h, 0)*Xp^3;
assert Evaluate(f5, Xp) eq Yp^2;

fI, L := IntegralModelPolynomial(f5);
C := HyperellipticCurve(fI);
J := Jacobian(C);

D := J![X-Xp, L*Yp];
assert Order(D) eq 12;
Dbuilt := 4*D;
assert Order(Dbuilt) eq 3;

G, phi := TorsionSubgroup(J);
invs := Invariants(G);

extra_found := false;
Pextra := J!0;
gextra := G!0;
for g in G do
    P := phi(g);
    if Order(P) eq 3 and P ne Dbuilt and P ne -Dbuilt then
        extra_found := true;
        Pextra := P;
        gextra := g;
        break;
    end if;
end for;

assert extra_found;
assert 3*Pextra eq J!0;
assert Pextra ne J!0;
assert Pextra ne Dbuilt;
assert Pextra ne -Dbuilt;

print "M(2,12)+extra independent 3 hit";
print "z", z, "r", r, "a", a, "w", w;
print "odd_integral_model", fI;
print "torsion", invs;
print "D_order", Order(D);
print "built_in_3 = 4D =", Dbuilt;
print "extra_3 =", Pextra;
print "abstract_extra_generator", gextra;
print "checks", "3*extra", 3*Pextra, "extra_not_builtin", Pextra ne Dbuilt and Pextra ne -Dbuilt;

quit;
