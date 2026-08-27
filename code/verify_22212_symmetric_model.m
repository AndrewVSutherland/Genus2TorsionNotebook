//////////////////////////////////////////////////////////////////////
// Verify the symmetric model for the (2,2,2,12) construction.
//
// The script has two deliberately small parts:
//   1. exact polynomial identities modulo the two surface equations;
//   2. exact torsion computations for the three known rational points.
//
// It does not construct a function-field Jacobian or run a generic
// torsion computation, so its memory use stays modest.
//
// Run:
//   magma -b code/verify_22212_symmetric_model.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

Q := Rationals();
R<x,y,z,u,v,T> := PolynomialRing(Q,6);

surface2 := x^2+y^2+z^2-u^2-v^2;
surface4 := x^4+y^4+z^4-u^4-v^4;
I := ideal<R | surface2, surface4>;

p := x^2;
q := y^2;
r := z^2;
b := u^2;
c := v^2;
d := p+q-b;
e := d+q;

// Elementary symmetric identities.
assert p+q+r-b-c in I;
assert p*q+p*r+q*r-b*c in I;
assert d-(c-r) in I;

// The two sextics in the original A,B,C square cover, on
// (A,B,C)=(x,u,y).
F_literal :=
    p^2*b-2*p*b^2+b^3-p^2*q+3*p*b*q-2*b^2*q-p*q^2+b*q^2;
G_literal :=
    p^3-2*p^2*b+p*b^2+2*p^2*q-3*p*b*q+b^2*q
    +2*p*q^2-2*b*q^2+q^3;
F_factored := d*(p-b)*(b-q);
G_factored := F_factored+d^3;
assert F_literal eq F_factored;
assert G_literal eq G_factored;
assert F_factored-d^2*r in I;
assert G_factored-d^2*c in I;

// Identities used by the order-four halving test.
assert (b-p)*(c-p)-q*r in I;
assert (b-q)*(c-q)-p*r in I;
assert (b-r)*(c-r)-p*q in I;

// Constant difference of the two cubic factors: this gives order 3.
Gcubic := (T-p)*(T-q)*(T-r);
Hcubic := T*(T-b)*(T-c);
assert Gcubic-Hcubic+p*q*r in I;

// Projective (s,m,n) and exact factorization of the five old branch
// parameters.
S0 := p*e;
M0 := 2*(q-b)*e;
N0 := 2*p*d;
lambda := 2*p*e;

// In particular these are exactly the three rational formulas in the
// question after clearing their common denominators.
M_literal := -2*p*b+2*b^2+2*p*q-6*b*q+4*q^2;
assert M_literal eq M0;

B_old := [
    2*S0^2-S0*N0,
    2*S0^2+S0*M0-2*S0*N0-M0*N0,
    2*S0^2+S0*M0-S0*N0-M0*N0,
    -M0*N0,
    4*S0^2-4*S0*N0-M0*N0
];
B_factored := [
    lambda*p*q,
    lambda*(b-p)*d,
    lambda*b*d,
    2*lambda*(b-q)*d,
    2*lambda*(p*q-d^2)
];
assert &and[B_old[i]-B_factored[i] in I : i in [1..5]];

// After removing the common square and putting R0=d*(T-b), the old
// sextic is d^6 times the symmetric sextic.  This also certifies that
// there is no hidden quadratic twist.
R0 := d*(T-b);
old_clean :=
    R0*(R0+p*q)*(R0+(b-p)*d)*(R0+b*d)
      *(R0+(b-q)*d)*(R0+p*q-d^2);
symmetric :=
    T*(T-p)*(T-q)*(T-r)*(T-b)*(T-c);
assert old_clean-d^6*symmetric in I;

printf "SYMBOLIC_IDENTITIES_VERIFIED\n";

//////////////////////////////////////////////////////////////////////
// Known rational specializations.
//////////////////////////////////////////////////////////////////////

Z := Integers();
P1<t> := PolynomialRing(Q);

function PrimitiveTriple(vals)
    ints := [Z!a : a in vals];
    gg := GCD([Abs(a) : a in ints]);
    ints := [a div gg : a in ints];
    if ints[1] lt 0 then
        ints := [-a : a in ints];
    end if;
    return ints;
end function;

function CurveFromSMN(ss,mm,nn)
    avec := [1,1,1,2,2];
    bvec := [
        2*ss^2-ss*nn,
        2*ss^2+ss*mm-2*ss*nn-mm*nn,
        2*ss^2+ss*mm-ss*nn-mm*nn,
        -mm*nn,
        4*ss^2-4*ss*nn-mm*nn
    ];
    ff := &*[avec[i]+bvec[i]*t : i in [1..5]];
    assert Degree(ff) eq 5 and Discriminant(ff) ne 0;
    return HyperellipticCurve(ff);
end function;

points := [
    <[408,143,1015,437,1013], [336396,-689185,-166464]>,
    <[120,143,266,218,241], [2208,-8303,-7200]>,
    <[16660,78793,21456,78644,27593],
       [1500518600,253638081,138777800]>
];

for item in points do
    coords := item[1];
    expected_smn := item[2];
    xx := coords[1];
    yy := coords[2];
    zz := coords[3];
    uu := coords[4];
    vv := coords[5];

    assert xx^2+yy^2+zz^2 eq uu^2+vv^2;
    assert xx^4+yy^4+zz^4 eq uu^4+vv^4;

    dd := xx^2+yy^2-uu^2;
    ee := dd+yy^2;
    smn := PrimitiveTriple([
        xx^2*ee,
        2*(yy^2-uu^2)*ee,
        2*xx^2*dd
    ]);
    assert smn eq expected_smn;

    fsym :=
        t*(t-xx^2)*(t-yy^2)*(t-zz^2)*(t-uu^2)*(t-vv^2);
    assert Degree(fsym) eq 6 and Discriminant(fsym) ne 0;
    Csym := HyperellipticCurve(fsym);
    Csmall := SimplifiedModel(Csym);

    Tor := TorsionSubgroup(Jacobian(Csmall));
    invs := Invariants(Tor);
    assert invs eq [2,2,2,12] and #Tor eq 96;

    Cold := CurveFromSMN(smn[1],smn[2],smn[3]);
    assert G2Invariants(Csym) eq G2Invariants(Cold);

    printf "POINT %o  SMN %o  TORSION %o  SAME_G2 true\n",
           coords,smn,invs;
end for;

printf "ALL_THREE_SYMMETRIC_MODELS_VERIFIED\n";
quit;
