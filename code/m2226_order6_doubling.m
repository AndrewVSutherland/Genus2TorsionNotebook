//////////////////////////////////////////////////////////////////////
//  M(2,2,2,6): standard quintic model and tests for order-6 classes
//  lying on the theta-doubling locus.
//
//  Start from the cleaner model in Theorem 9.4:
//
//      y^2 = x (x + 2s^2 - sn)
//              (x + 2s^2 + sm - 2sn - mn)
//              (x + 2s^2 + sm - sn - mn)
//              (2x - mn)
//              (2x + 4s^2 - 4sn - mn).
//
//  Sending the Weierstrass point x=0 to infinity by
//
//      X = 1/x,  Y = y/x^3
//
//  gives an odd-degree model
//
//      Y^2 = prod_i L_i(X),       P = (0,2),
//
//  where infinity is the old Weierstrass point x=0 and P is one of the
//  old points at infinity.  The class P-infinity has order 6.
//
//  The six order-6 classes represented by P plus a Weierstrass point are
//  not considered here.  The other ten are
//
//      (P-infinity) + (W_i + W_j - 2 infinity),  1 <= i < j <= 5.
//
//  This script writes their Mumford representatives in a compact form
//  and computes the discriminant equations for the condition that the
//  representative divisor is 2Q for a point Q on C.
//
//  Typical run from torsion_jac:
//
//      magma code/m2226_order6_doubling.m
//////////////////////////////////////////////////////////////////////

Q := Rationals();
R<s,m,n> := PolynomialRing(Q, 3);
PX<X> := PolynomialRing(R);
K := FieldOfFractions(R);
PK<XK> := PolynomialRing(K);

// Finite Weierstrass factors after X=1/x.  The constants are separated
// so that the conic through P and two Weierstrass points has a simple
// expression even for the two factors with leading coefficient 2.
A := [ R!1, R!1, R!1, R!2, R!2 ];
B := [
    2*s^2 - s*n,
    2*s^2 + s*m - 2*s*n - m*n,
    2*s^2 + s*m - s*n - m*n,
    -m*n,
    4*s^2 - 4*s*n - m*n
];
L := [ PX!A[i] + PX!B[i]*X : i in [1..5] ];
f := &*L;
assert Coefficient(f, 0) eq 4;

function ProductExcept(seq, omitted)
    out := Universe(seq)!1;
    for i in [1..#seq] do
        if i notin omitted then
            out *:= seq[i];
        end if;
    end for;
    return out;
end function;

function PrimitivePolynomial(g)
    if g eq 0 then
        return g;
    end if;
    coeffs := Coefficients(g);
    den_lcm := LCM([ Denominator(c) : c in coeffs ]);
    h := den_lcm*g;
    nums := [ Integers()!c : c in Coefficients(h) ];
    cont := GCD([ Abs(c) : c in nums | c ne 0 ]);
    if cont ne 0 then
        h /:= cont;
    end if;
    if LeadingCoefficient(h) lt 0 then
        h := -h;
    end if;
    return h;
end function;

function PairRawU(i, j)
    // Let beta_i be the root of L_i.  The quadratic through
    // P=(0,2), W_i=(beta_i,0), and W_j=(beta_j,0) is
    //
    //     ell_ij = 2 L_i L_j/(A_i A_j).
    //
    // If G is the product of the other three L's, then the two residual
    // intersection points have x-coordinates cut out by
    //
    //     B_ij(X) = (G - 4 L_i L_j/(A_i^2 A_j^2))/X.
    //
    // The Mumford U-polynomial is the monic normalization of B_ij.
    G := ProductExcept(L, {i,j});
    numerator := G - (Q!4)/(A[i]^2*A[j]^2)*L[i]*L[j];
    ok, Bij := IsDivisibleBy(numerator, X);
    assert ok;
    assert Degree(Bij) eq 2;
    return Bij;
end function;

function PairEll(i, j)
    return (Q!2)/(A[i]*A[j])*L[i]*L[j];
end function;

function ToPK(poly)
    out := PK!0;
    for k in [0..Degree(poly)] do
        out +:= K!Coefficient(poly, k)*XK^k;
    end for;
    return out;
end function;

function PairMonicU(i, j)
    Bij := ToPK(PairRawU(i, j));
    return Bij/LeadingCoefficient(Bij);
end function;

function PairV(i, j)
    Uij := PairMonicU(i, j);
    ell := ToPK(PairEll(i, j));
    _, Vij := Quotrem(-ell, Uij);
    return Vij;
end function;

function PairDiscriminant(i, j)
    return PrimitivePolynomial(Discriminant(PairRawU(i, j)));
end function;

procedure PrintQuinticModel()
    print "Odd-degree model after X=1/x, Y=y/x^3:";
    print "Y^2 =", f;
    print "P = (0, 2), infinity is the old x=0 Weierstrass point.";
    print "";
    for i in [1..5] do
        print "L", i, "=", L[i];
    end for;
    print "";
end procedure;

procedure PrintPairData()
    print "For each pair i<j:";
    print "  ell_ij = 2*L_i*L_j/(A_i*A_j)";
    print "  raw U_ij = (prod_{k notin {i,j}} L_k - 4*L_i*L_j/(A_i^2*A_j^2))/X";
    print "  Mumford representative is [monic(raw U_ij), -ell_ij mod monic(raw U_ij)].";
    print "  The theta-doubling condition is Disc_X(raw U_ij)=0.";
    print "";

    for i in [1..4] do
        for j in [i+1..5] do
            print "PAIR", [i,j];
            print "raw_U =", PairRawU(i, j);
            print "ell =", PairEll(i, j);
            print "V is the remainder of -ell modulo monic(raw_U).";
            disc := PairDiscriminant(i, j);
            print "disc =", disc;
            if disc eq 0 then
                print "factorization = identically zero";
            else
                print "factorization =", Factorization(disc);
            end if;
            print "";
        end for;
    end for;
end procedure;

function EvalRationalPolynomial(g, vals)
    return Evaluate(g, vals);
end function;

function EvalPX(poly, vals, xvar)
    Qx := Parent(xvar);
    out := Qx!0;
    for k in [0..Degree(poly)] do
        out +:= Q!EvalRationalPolynomial(Coefficient(poly, k), vals)*xvar^k;
    end for;
    return out;
end function;

procedure VerifySpecialization(vals)
    Qx<x> := PolynomialRing(Q);
    Lq := [ EvalPX(Li, vals, x) : Li in L ];
    fq := &*Lq;

    if Discriminant(fq) eq 0 then
        print "Skipping singular specialization", vals;
        return;
    end if;

    C := HyperellipticCurve(fq);
    J := Jacobian(C);
    g := J![x, Q!2];
    assert 6*g eq J!0;
    assert 3*g ne J!0;

    Aq := [ Q!Evaluate(Ai, vals) : Ai in A ];
    Bq := [ Q!Evaluate(Bi, vals) : Bi in B ];

    for i in [1..4] do
        for j in [i+1..5] do
            if Bq[i] eq 0 or Bq[j] eq 0 then
                print "Skipping pair with Weierstrass point at X=infinity after specialization", vals, [i,j];
                continue;
            end if;

            beta_i := -Aq[i]/Bq[i];
            beta_j := -Aq[j]/Bq[j];
            Tij := J![(x-beta_i)*(x-beta_j), Q!0];

            Uq_raw := EvalPX(PairRawU(i, j), vals, x);
            Uq := Uq_raw/LeadingCoefficient(Uq_raw);
            ellq := EvalPX(PairEll(i, j), vals, x);
            _, Vq := Quotrem(-ellq, Uq);
            Dij := J![Uq, Vq];

            if Dij ne g + Tij then
                print "FAILED specialization", vals, "pair", [i,j];
                print "fq =", fq;
                print "Uq =", Uq;
                print "Vq =", Vq;
                print "Dij =", Dij;
                print "g+Tij =", g+Tij;
                error "Pair formula check failed.";
            end if;

            ord := Order(g + Tij);
            if i eq 3 and j eq 5 then
                assert ord eq 3;
            else
                assert ord eq 6;
            end if;
        end for;
    end for;

    print "Verified specialization", vals;
end procedure;

PrintQuinticModel();
PrintPairData();

VerifySpecialization([ Q!1, Q!3, Q!5 ]);
VerifySpecialization([ Q!2, Q!-1, Q!3 ]);
