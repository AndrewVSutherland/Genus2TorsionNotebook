//////////////////////////////////////////////////////////////////////
//  Derive the automorphism/split condition on the eps=+1 extra-root branch.
//
//  On eps=+1, a rational extra root r is a root of
//
//      Qp = (b+3)*x^2 + (a-3)*x + 2.
//
//  Let s be the second root of Qp.  The Mobius involution swapping
//
//      0 <-> r,    s <-> infinity
//
//  is tau(x)=s*(r-x)/(s-x).  The genus-2 curve has an extra involution if
//  tau preserves the residual quadratic
//
//      Qm = 2*x^2 + (b-3)*x + (a+3).
//
//  This is a concrete Humbert/decomposable condition: it gives degree-2
//  elliptic subcovers and hence a split Jacobian.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
R<b,r> := PolynomialRing(Q, 2);
K := FieldOfFractions(R);
P<x> := PolynomialRing(K);

a := 3 - (b+3)*r - 2/r;
s := 2/((b+3)*r);
Qm := 2*x^2 + (b-3)*x + (a+3);
tau := s*(r-x)/(s-x);

T := P!Numerator((s-x)^2 * Evaluate(Qm, tau));
// T and Qm have the same leading/constant scaling iff the cross-products
// of coefficients vanish.
cT := [Coefficient(T,i) : i in [0..2]];
cQ := [Coefficient(Qm,i) : i in [0..2]];
eqs := [];
for i in [1..3] do
    for j in [i+1..3] do
        Append(~eqs, Numerator(cT[i]*cQ[j] - cT[j]*cQ[i]));
    end for;
end for;
eqsR := [R!Numerator(e) : e in eqs];
g := GCD(eqsR);
gR := R!Numerator(g);
fac := Factorization(gR);


// Alternate pairing: 0 <-> s and r <-> infinity.
tau2 := r*(s-x)/(r-x);
T2 := P!Numerator((r-x)^2 * Evaluate(Qm, tau2));
cT2 := [Coefficient(T2,i) : i in [0..2]];
eqs2 := [];
for i in [1..3] do
    for j in [i+1..3] do
        Append(~eqs2, Numerator(cT2[i]*cQ[j] - cT2[j]*cQ[i]));
    end for;
end for;
eqs2R := [R!Numerator(e) : e in eqs2];

function EvalList(polys, pt)
    return [Evaluate(h, pt) : h in polys];
end function;

print "eps=+1 extra-root automorphism condition";
print "a = 3 - (b+3)*r - 2/r";
print "second Qp root s = 2/((b+3)*r)";
print "tau(x)=s*(r-x)/(s-x)";
print "pairing A equations =", eqsR;
print "pairing B equations =", eqs2R;
print "gcd necessary factor =", gR;
print "factorization", fac;
print "known hit r=21,b=-23/7 A", EvalList(eqsR, <Q!(-23/7),Q!21>), "B", EvalList(eqs2R, <Q!(-23/7),Q!21>);
print "known hit r=-1/3,b=-23/7 A", EvalList(eqsR, <Q!(-23/7),Q!(-1/3)>), "B", EvalList(eqs2R, <Q!(-23/7),Q!(-1/3)>);
print "known hit r=4/3,b=-13/7 A", EvalList(eqsR, <Q!(-13/7),Q!(4/3)>), "B", EvalList(eqs2R, <Q!(-13/7),Q!(4/3)>);
print "generic simple probe r=2,b=0 A", EvalList(eqsR, <Q!0,Q!2>), "B", EvalList(eqs2R, <Q!0,Q!2>);

// Reduce modulo 19 and compare the finite red/irreducible samples in eps=1,r=2.
F := GF(19);
R19<bb,rr> := PolynomialRing(F, 2);
phi := hom<R -> R19 | bb,rr>;
eqs19 := [phi(h) : h in eqsR];
eqs219 := [phi(h) : h in eqs2R];
print "mod19 equations at r=2,b=13 A", [Evaluate(h, <F!13,F!2>) : h in eqs19], "B", [Evaluate(h, <F!13,F!2>) : h in eqs219];
print "mod19 equations at r=2,b=3 A", [Evaluate(h, <F!3,F!2>) : h in eqs19], "B", [Evaluate(h, <F!3,F!2>) : h in eqs219];

quit;
