//////////////////////////////////////////////////////////////////////
//  Lightweight p=19 component diagnostic for eps=1, r=2.
//
//  Avoids primary decomposition.  It saturates the fixed chart, tests
//  primality when possible, enumerates F_19 points, and compares local
//  tangent dimensions at reducible and irreducible Frobenius points.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

F := GF(19);
Z := Integers();
R<b,L,U,v> := PolynomialRing(F, 4);
PR<x> := PolynomialRing(R);

eps := F!1;
r := F!2;
a := (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;
h6 := 1 + a*x + b*x^2 + x^3;
f := h6^2 - (x-1)^6;
c1 := Coefficient(f, 1);
c2 := Coefficient(f, 2);
c3 := Coefficient(f, 3);
c4 := Coefficient(f, 4);
c5 := Coefficient(f, 5);
B := c5*L^2 + 3*U;
Delta := 4*c4*L^2 + 12*(U^2 + v^2) - B^2;
F3 := B*Delta + 16*v^3 - 8*c3*L^2 - 8*U^3 - 48*U*v^2;
F2 := Delta^2 + 64*B*v^3 - 64*c2*L^2 - 192*(U^2*v^2 + v^4);
F1 := Delta*v^3 - 4*c1*L^2 - 12*U*v^4;
q := x^2 + U*x + v^2;

boundary := L*v*(U^2 - 4*v^2)*(b+3)*Evaluate(h6, F!1)*Discriminant(f)*Discriminant(q);
I := ideal<R | F1,F2,F3>;
Isat := Saturation(I, ideal<R | boundary>);

function EvalPoly(poly, pt)
    return Evaluate(poly, <pt[1],pt[2],pt[3],pt[4]>);
end function;

function SpecializePoly(poly, pt)
    PF<X> := PolynomialRing(F);
    coeffs := [EvalPoly(Coefficient(poly,i), pt) : i in [0..Degree(poly)]];
    return PF!coeffs;
end function;

function GoodPoint(pt)
    aa := EvalPoly(a, pt);
    ff := SpecializePoly(f, pt);
    hh := SpecializePoly(h6, pt);
    qq := SpecializePoly(q, pt);
    if Degree(ff) ne 5 or Discriminant(ff) eq 0 then
        return false, false, [], ff, Parent(ff)!0;
    end if;
    if Evaluate(hh, r) eq 0 or Evaluate(ff, r) ne 0 then
        return false, false, [], ff, Parent(ff)!0;
    end if;
    if Discriminant(qq) eq 0 or Degree(GCD(qq, ff)) gt 0 then
        return false, false, [], ff, Parent(ff)!0;
    end if;
    C := HyperellipticCurve(ff);
    J := Jacobian(C);
    A, phi := AbelianGroup(J);
    Lp := LPolynomial(C);
    fac := Factorization(Lp);
    irred := #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4;
    return true, irred, Invariants(A), ff, Lp;
end function;

function PointInIdeal(J, pt)
    return &and[EvalPoly(g, pt) eq 0 : g in Basis(J)];
end function;

function TangentDim(J, pt)
    rows := [];
    for g in Basis(J) do
        Append(~rows, [EvalPoly(Derivative(g,i), pt) : i in [1..4]]);
    end for;
    M := Matrix(F, rows);
    return 4 - Rank(M);
end function;

known_red := <F!13,F!6,F!0,F!8>;
known_irred := <F!3,F!8,F!3,F!14>;

print "p=19 light component diagnostic for eps=1, r=2";
dim, degs := Dimension(Isat);
print "saturated_dimension", dim, "component_degrees", degs;
try
    print "degree", Degree(Isat);
catch e
    print "degree unavailable", e`Object;
end try;
try
    print "is_prime", IsPrime(Isat);
catch e
    print "is_prime unavailable", e`Object;
end try;
print "known_red_in_Isat", PointInIdeal(Isat, known_red),
      "tangent_dim", TangentDim(Isat, known_red);
print "known_irred_in_Isat", PointInIdeal(Isat, known_irred),
      "tangent_dim", TangentDim(Isat, known_irred);

points := [];
for bb in F do
    for LL in F do
        for UU in F do
            for vv in F do
                pt := <bb,LL,UU,vv>;
                if PointInIdeal(Isat, pt) then
                    Append(~points, pt);
                end if;
            end for;
        end for;
    end for;
end for;

good := 0;
irred := 0;
red := 0;
by_b := AssociativeArray();
sample_irred := [];
sample_red := [];
for pt in points do
    ok, isirred, invs, ff, Lp := GoodPoint(pt);
    if not ok then
        continue;
    end if;
    good +:= 1;
    key := Z!pt[1];
    if not IsDefined(by_b, key) then
        by_b[key] := <0,0,0>;
    end if;
    rec := by_b[key];
    if isirred then
        irred +:= 1;
        by_b[key] := <rec[1]+1, rec[2]+1, rec[3]>;
        if #sample_irred lt 5 then
            Append(~sample_irred, <pt,invs,Lp>);
        end if;
    else
        red +:= 1;
        by_b[key] := <rec[1]+1, rec[2], rec[3]+1>;
        if #sample_red lt 5 then
            Append(~sample_red, <pt,invs,Lp>);
        end if;
    end if;
end for;

print "F19_points_on_Isat", #points;
print "good", good, "irred_L", irred, "red_L", red;
print "by_b <good,irred,red>", by_b;
print "sample_irred", sample_irred;
print "sample_red", sample_red;

quit;
