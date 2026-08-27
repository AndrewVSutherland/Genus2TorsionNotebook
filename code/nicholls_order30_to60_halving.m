//////////////////////////////////////////////////////////////////////
//  Audit Nicholls, Example 3.1.8, as a possible order-60 seed.
//
//  The script verifies:
//    * the exact rational torsion subgroup is Z/30Z;
//    * the sextic has irreducible factor degrees 2+4, giving one
//      nonzero rational 2-class T2;
//    * T2 is 15 times a generator of the order-30 subgroup;
//    * at the good prime 5, J(F_5) is cyclic of order 30 and T2 is
//      not a double, so no rational half of T2 can exist;
//    * Magma's exact IsDivisibleBy agrees with the local obstruction;
//    * the good reduction at 13 gives the D4/root-power geometric
//      simplicity certificate used elsewhere in this repository.
//
//  Run from the repository root with
//
//    magma code/nicholls_order30_to60_halving.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

f := x^6 - (Q!16/3)*x^5 + (Q!70/9)*x^4
     + (Q!131/27)*x^2 + (Q!16/27)*x + Q!64/81;

// The square change Y=9y gives an integral model of the same curve.
fI := P!(81*f);
assert fI eq 81*x^6 - 432*x^5 + 630*x^4
               + 393*x^2 + 48*x + 64;

q2 := x^2 + Q!1/3;
q4 := x^4 - (Q!16/3)*x^3 + (Q!67/9)*x^2
      + (Q!16/9)*x + Q!64/27;
assert fI eq 81*q2*q4;
assert IsIrreducible(q2) and IsIrreducible(q4);

fac := Factorization(fI);
degrees := [Degree(pe[1]) : pe in fac];
Sort(~degrees);
assert degrees eq [2,4];

C := HyperellipticCurve(fI);
J := Jacobian(C);

// Exact rational torsion.
G, phi := TorsionSubgroup(J);
invs := Invariants(G);
assert invs eq [30];
P30 := phi(G.1);
assert Order(P30) eq 30;

// The unique nonzero rational 2-class furnished by the 2+4 factorization.
T2 := J![q2, Q!0];
assert Order(T2) eq 2;
assert 15*P30 eq T2;

//////////////////////////////////////////////////////////////////////
//  Local halving obstruction at p=5.
//////////////////////////////////////////////////////////////////////

p := 5;
F := GF(p);
PF<z> := PolynomialRing(F);
f5 := PF![F!Coefficient(fI,i) : i in [0..Degree(fI)]];
q2_5 := PF![F!Coefficient(q2,i) : i in [0..Degree(q2)]];
assert Discriminant(f5) ne 0;

C5 := HyperellipticCurve(f5);
J5 := Jacobian(C5);
A5, amap5 := AbelianGroup(J5);
assert Invariants(A5) eq [30];

T2_5 := J5![q2_5, F!0];
assert Order(T2_5) eq 2;
t2_abstract := T2_5 @@ amap5;
doubles5 := sub<A5 | [2*A5.i : i in [1..Ngens(A5)]]>;
local_half_exists := t2_abstract in doubles5;
assert not local_half_exists;

// Independent exact cross-check over Q; the local obstruction already proves this.
global_half_exists, half := IsDivisibleBy(T2, 2);
assert not global_half_exists;

//////////////////////////////////////////////////////////////////////
//  Geometric-simplicity certificate at p=13.
//////////////////////////////////////////////////////////////////////

PT<T> := PolynomialRing(Q);

function FrobeniusPolynomial(C0, ell)
    ef := EulerFactor(C0, ell);
    d := Degree(ef);
    return &+[Q!Coefficient(ef,i)*T^(d-i) : i in [0..d]];
end function;

Phi13 := FrobeniusPolynomial(C, 13);
assert Phi13 eq T^4 - 4*T^3 + 6*T^2 - 52*T + 169;
fac13 := Factorization(Phi13);
assert #fac13 eq 1 and fac13[1][2] eq 1
       and Degree(fac13[1][1]) eq 4;

Gal13 := GaloisGroup(Phi13);
gal_desc := TransitiveGroupDescription(Gal13);
assert Order(Gal13) eq 8 and gal_desc eq "D(4)";

K13<pi13> := NumberField(Phi13);
power_degrees := [Degree(MinimalPolynomial(pi13^n)) : n in [2..12]];
assert &and[d eq 4 : d in power_degrees];

print "Nicholls Example 3.1.8 order-30 seed audit";
print "rational model f =", f;
print "integral square-scaled model 81*f =", fI;
print "irreducible factor degrees =", degrees;
print "unique rational 2-class T2 =", T2;
print "torsion invariants =", invs;
print "15*P30 = T2 =", 15*P30 eq T2;
print "p=5 finite Jacobian invariants =", Invariants(A5);
print "T2 is a double in J(F_5) =", local_half_exists;
print "T2 is a double in J(Q) =", global_half_exists;
print "p=13 Frobenius polynomial =", Phi13;
print "p=13 Galois group order/description =", Order(Gal13), gal_desc;
print "degrees of minpoly(pi^n), n=2..12 =", power_degrees;
print "RESULT: exact torsion Z/30; p=5 obstructs halving; geometrically simple.";

quit;
