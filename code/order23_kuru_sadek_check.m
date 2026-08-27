//////////////////////////////////////////////////////////////////////
// Kuru-Sadek genus-2 order-23 family.
// Source: Quadratic torsion orders on Jacobian varieties, arXiv:2410.14455.
// The alpha formula uses beta^((g+1)/2); this reproduces their t=2 example.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
P<x> := PolynomialRing(Q);

function KuruSadek23(t)
    t := Q!t;
    beta := (t^2 + 1)^2/(4*t^2);
    sbeta := (t^2 + 1)/(2*t);
    s := (t^2 - 1)/(2*t);
    alpha := beta - s^5/(beta*sbeta);
    lambda := (alpha - 1)^4/((alpha - beta)^2*alpha);
    num := x^3*(x - alpha)^2 - (x - 1)*((x - 1)^4 - lambda*(x - beta)^2*x);
    den := 2*(x - alpha)*(x - beta);
    A, rem := Quotrem(num, den);
    assert rem eq 0;
    return A^2 - lambda*x^4*(x - 1);
end function;

function PrimitiveIntegralSquareModel(f)
    den := LCM([ Denominator(Coefficient(f, i)) : i in [0..Degree(f)] ]);
    F := den*f;
    cont := GCD([ Integers()!Coefficient(F, i) : i in [0..Degree(F)] ]);
    scale := den/cont;
    ok, root := IsSquare(Numerator(scale)*Denominator(scale));
    assert ok;
    return F/cont, scale;
end function;

f := KuruSadek23(2);
F, scale := PrimitiveIntegralSquareModel(f);
print "f", f;
print "square scale to primitive model", scale;
print "primitive F", F;
print "Discriminant zero", Discriminant(F) eq 0;
C := HyperellipticCurve(F);
J := Jacobian(C);
for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
    fp := ChangeRing(F, GF(p));
    if Discriminant(fp) eq 0 then continue; end if;
    Lp := LPolynomial(ChangeRing(C, GF(p)));
    fac := Factorization(Lp);
    if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
        print "irreducible Lp", p, Lp;
        break;
    end if;
end for;
T, mp := TorsionSubgroup(J);
print "TorsionSubgroup", Invariants(T);
