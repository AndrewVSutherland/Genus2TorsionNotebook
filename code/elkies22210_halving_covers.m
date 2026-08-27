//////////////////////////////////////////////////////////////////////
//  Exact halving covers over Elkies's Clebsch--Klein surface.
//
//  The source curve is
//
//      C_r : y^2 = x * Product_{i=1}^5 (x-r_i^2),
//
//  with Sum r_i = Sum r_i^3 = 0.  On the smooth open, the six
//  finite branch roots are 0,r_1^2,...,r_5^2 and J[2](Q) is full.
//  There are two S_5-orbits of nonzero 2-torsion classes:
//
//      T_01 = {0,r_1^2},          T_12 = {r_1^2,r_2^2}.
//
//  This file derives and certifies the degree-16 rational-halving
//  covers for representatives of both orbits.  It also checks the
//  criterion against Magma's exact IsDivisibleBy on all 15 classes of
//  Elkies's printed source and on a positive control.
//
//  Run from torsion_jac with
//
//      magma -b code/elkies22210_halving_covers.m
//////////////////////////////////////////////////////////////////////

Q := Rationals();
Z := Integers();

R<r1,r2,r3,r4,r5> := PolynomialRing(Q, 5);
rs := [r1,r2,r3,r4,r5];
aa := [ri^2 : ri in rs];

ck_linear := &+rs;
ck_cubic := &+[ri^3 : ri in rs];
boundary := &*rs * &*[ aa[i]-aa[j] : i,j in [1..5] | i lt j ];

// Orbit {0,r_1^2}.  After removing square factors from the exact
// Stoll/Zarhin radicands, the cover is simply
//
//      z_j^2 = r_1^2-r_j^2,  j=2,3,4,5.
orbit01 := [ aa[1]-aa[j] : j in [2..5] ];

// Orbit {r_1^2,r_2^2}.  Again remove the visible square factors r_2^2
// from the first exact radicand and r_1^2 from the other three.
orbit12 := [
    -(aa[1]-aa[3])*(aa[1]-aa[4])*(aa[1]-aa[5]),
    (aa[3]-aa[2])*(aa[1]-aa[4])*(aa[1]-aa[5]),
    (aa[4]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[5]),
    (aa[5]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[4])
];

K := FieldOfFractions(R);
KT<T> := PolynomialRing(K);

// If A and B are two branch roots and c runs over the other four,
// put
//
//   t=(x-B)/(x-A), D=Product_c(A-c), lambda_c=(c-B)/(c-A).
//
// The transformed monic odd model has roots 0 and D*lambda_c.
// Stoll's Kummer map therefore says that B-A is twice a rational
// divisor precisely when every -D*lambda_c is a square.
function ExactRadicands(A, B, remaining)
    A := K!A;
    B := K!B;
    remaining := [K!c : c in remaining];
    D := &*[A-c : c in remaining];
    vals := [ -D*(c-B)/(c-A) : c in remaining ];
    return D, vals;
end function;

// Certify the fractional-linear transformation identity.  If the
// original even model is monic, then after
//
//   X=D*(x-B)/(x-A),
//   Y=D^2*(A-B)^2*y/(x-A)^3,
//
// it is Y^2=X*Product_c(X-D*lambda_c).
procedure CheckTransformIdentity(A, B, remaining)
    A := K!A;
    B := K!B;
    remaining := [K!c : c in remaining];
    D := &*[A-c : c in remaining];
    lambdas := [(c-B)/(c-A) : c in remaining];
    x_of_t := (T*A-B)/(T-1);
    lhs := (T-1)^6 * (x_of_t-A)*(x_of_t-B)
        * &*[x_of_t-c : c in remaining];
    rhs := T*(A-B)^2*D*&*[T-lam : lam in lambdas];
    assert lhs eq rhs;

    // Scaling X=D*T makes the odd polynomial monic.
    PX<X> := PolynomialRing(K);
    lhs_monic := X*&*[X-D*lam : lam in lambdas];
    rhs_monic := D^5 * ( (X/D) * &*[X/D-lam : lam in lambdas] );
    assert lhs_monic eq rhs_monic;
end procedure;

// Symbolic derivation for T_01.
D01, exact01 := ExactRadicands(0, aa[1], aa[2..5]);
prod_r_2_5 := &*rs[2..5];
assert D01 eq K!(prod_r_2_5^2);
for j in [2..5] do
    // -D*lambda_j = (Product_{k=2}^5 r_k / r_j)^2
    //                  * (r_1^2-r_j^2).
    assert exact01[j-1] eq
        K!((prod_r_2_5 div rs[j])^2 * orbit01[j-1]);
end for;
CheckTransformIdentity(0, aa[1], aa[2..5]);

// Symbolic derivation for T_12.
remaining12 := [R!0, aa[3], aa[4], aa[5]];
D12, exact12 := ExactRadicands(aa[1], aa[2], remaining12);
assert exact12[1] eq K!(aa[2]*orbit12[1]);
for j in [2..4] do
    assert exact12[j] eq K!(aa[1]*orbit12[j]);
end for;
for j in [3..5] do
    assert K!(orbit12[j-1]/orbit12[1]) eq
        K!(-(aa[j]-aa[2])/(aa[1]-aa[j]));
end for;
CheckTransformIdentity(aa[1], aa[2], remaining12);

// Rational two-parameter chart of the full labelled smooth CK open.
// Pair the coordinates as r1=x+y,r2=x-y,r3=z+w,r4=z-w and scale
// x=1, z=t.  A line through (-1,-t-2) on the resulting conic gives R.
CP<tc,mc> := PolynomialRing(Q,2);
ck_param := [
    1+tc*(tc+2)*mc,
    tc*mc*(mc-tc-2),
    -1+mc+tc*(tc+1)*mc^2,
    1+tc-mc-tc*mc^2,
    -(1+tc)*(1+tc*mc^2)
];
assert &+ck_param eq 0;
assert &+[u^3 : u in ck_param] eq 0;
xc := (ck_param[1]+ck_param[2])/2;
yc := (ck_param[1]-ck_param[2])/2;
zc := (ck_param[3]+ck_param[4])/2;
wc := (ck_param[3]-ck_param[4])/2;
assert zc eq tc*xc;
assert yc^2+tc*wc^2 eq
    xc^2*(1+tc)*(tc^2+3*tc+1);

// Conversely, on the conic recover m=(W+t+2)/(Y+1).  Substitution
// recovers [1+Y,1-Y,t+W,t-W,-2(1+t)] projectively.  The only excluded
// denominators are r1+r2=0 and r1=0, both on the CK boundary.
IP<ti,Yi,Wi> := PolynomialRing(Q,3);
FIP := FieldOfFractions(IP);
conic_inverse := Yi^2+ti*Wi^2-(1+ti)*(ti^2+3*ti+1);
mi := FIP!((Wi+ti+2)/(Yi+1));
inverse_param := [
    1+ti*(ti+2)*mi,
    ti*mi*(mi-ti-2),
    -1+mi+ti*(ti+1)*mi^2,
    1+ti-mi-ti*mi^2,
    -(1+ti)*(1+ti*mi^2)
];
normalized_r := [1+Yi,1-Yi,ti+Wi,ti-Wi,-2*(1+ti)];
inverse_scale := 1+ti*mi^2;
Iconic := ideal<IP | conic_inverse>;
for j in [1..5] do
    assert Numerator(2*inverse_param[j]
        -inverse_scale*normalized_r[j]) in Iconic;
end for;

// Direct unit-circle chart on orbit {0,r_1^2}.  Scale r_1=1 and
// parametrize each r_j^2+z_j^2=1.  On the linear CK equation, the
// cubic CK equation reduces to Sum r_j*z_j^2=0.
TP<t2,t3,t4,t5> := PolynomialRing(Q,4);
FP := FieldOfFractions(TP);
tpars := [FP!t2,FP!t3,FP!t4,FP!t5];
rpars := [(1-t^2)/(1+t^2) : t in tpars];
zpars := [2*t/(1+t^2) : t in tpars];
for j in [1..4] do
    assert rpars[j]^2+zpars[j]^2 eq 1;
end for;
chart_linear := 1+&+rpars;
chart_cubic := 1+&+[r^3 : r in rpars];
chart_weighted := &+[rpars[j]*zpars[j]^2 : j in [1..4]];
assert chart_cubic-chart_linear eq -chart_weighted;
chart_sum_recip := &+[1/(1+t^2) : t in tpars] - 3/2;
assert chart_linear eq 2*chart_sum_recip;
e2r := &+[rpars[i]*rpars[j] : i,j in [1..4] | i lt j];
e3r := &+[rpars[i]*rpars[j]*rpars[k] : i,j,k in [1..4] |
    i lt j and j lt k];
chart_den := &*[1+t^2 : t in tpars];
chart_product := &*[1-t^2 : t in tpars]-16;
assert chart_product eq -chart_den*(chart_linear+e2r+e3r);

// If s_j=t_j^2 and E_k=e_k(s_2,...,s_5), the two direct equations
// become simply E_4=E_1+5 and E_3=E_2-10.  This also gives an exact
// two-parameter pair-completion formula.  For pair sum/product A,B of
// s_2,s_3, the sum/product C,D of s_4,s_5 are the displayed functions.
spars := [t^2 : t in tpars];
E1s := &+spars;
E2s := &+[spars[i]*spars[j] : i,j in [1..4] | i lt j];
E3s := &+[spars[i]*spars[j]*spars[k] : i,j,k in [1..4] |
    i lt j and j lt k];
E4s := &*spars;
chart_symm_1 := E4s-E1s-5;
chart_symm_2 := E3s-E2s+10;
assert 2*chart_den*chart_sum_recip eq -3*chart_symm_1-chart_symm_2;
assert chart_product eq chart_symm_1-chart_symm_2;

A01 := spars[1]+spars[2];
B01 := spars[1]*spars[2];
C01 := spars[3]+spars[4];
D01pair := spars[3]*spars[4];
pair_den := (B01-1)*(B01-A01+1);
pair_C_num := B01^2-10*B01-A01^2-4*A01+5;
pair_D_num := (A01+6)*B01-(A01^2+5*A01+10);
assert TP!(pair_den*C01-pair_C_num) in ideal<TP |
    Numerator(chart_symm_1), Numerator(chart_symm_2)>;
assert TP!(pair_den*D01pair-pair_D_num) in ideal<TP |
    Numerator(chart_symm_1), Numerator(chart_symm_2)>;

// Zarhin's explicit half on the transformed monic odd model.
// For Y^2=X*Product_i(X+u_i^2), set e_i=e_i(u_1,...,u_4).
// The asserted congruence verifies that [q,alpha] is a divisor class;
// the subsequent exact Magma check verifies 2[q,alpha]=[(X),0].
S<u1,u2,u3,u4> := PolynomialRing(Q, 4);
PS<XX> := PolynomialRing(S);
us := [u1,u2,u3,u4];
e1 := &+us;
e2 := &+[us[i]*us[j] : i,j in [1..4] | i lt j];
e3 := &+[us[i]*us[j]*us[k] : i,j,k in [1..4] | i lt j and j lt k];
e4 := &*us;
qhalf := XX^2-e2*XX+e4;
alphahalf := (e1*e2-e3)*XX-e1*e4;
fzarhin := XX*&*[XX+ui^2 : ui in us];
assert (alphahalf^2-fzarhin) mod qhalf eq 0;

Pq<x> := PolynomialRing(Q);
uz := [Q!1,Q!2,Q!3,Q!4];
e1q := &+uz;
e2q := &+[uz[i]*uz[j] : i,j in [1..4] | i lt j];
e3q := &+[uz[i]*uz[j]*uz[k] : i,j,k in [1..4] | i lt j and j lt k];
e4q := &*uz;
fz := x*&*[x+ui^2 : ui in uz];
Cz := HyperellipticCurve(fz);
Jz := Jacobian(Cz);
Hz := Jz![x^2-e2q*x+e4q, (e1q*e2q-e3q)*x-e1q*e4q];
Tz := Jz![x,Q!0];
assert 2*Hz eq Tz;
assert Order(Hz) eq 4;

// Return the F_2-rank of rational squareclasses.  Rank m certifies
// that adjoining the m displayed square roots has degree 2^m.
function SquareclassRank(vals)
    vals := [Q!v : v in vals];
    primes := [];
    for v in vals do
        assert v ne 0;
        n := Numerator(v);
        d := Denominator(v);
        for pe in Factorization(Abs(n)) cat Factorization(d) do
            if pe[1] notin primes then
                Append(~primes, pe[1]);
            end if;
        end for;
    end for;
    Sort(~primes);
    F2 := GF(2);
    M := ZeroMatrix(F2, #vals, #primes+1);
    for i in [1..#vals] do
        v := vals[i];
        if v lt 0 then
            M[i,1] := 1;
        end if;
        for j in [1..#primes] do
            p := primes[j];
            M[i,j+1] := F2!(Valuation(Numerator(v),p)
                - Valuation(Denominator(v),p));
        end for;
    end for;
    return Rank(M), M, primes;
end function;

function EvalList(polys, vals)
    return [Q!Evaluate(f, vals) : f in polys];
end function;

// Elkies's simplest Clebsch--Klein point.  The four squareclasses in
// each orbit representative are independent, proving that each
// universal cover has generic degree 16 (specialization cannot have
// larger multiquadratic degree than the generic cover).
sample := [Q!1,Q!-8,Q!-7,Q!5,Q!9];
assert Evaluate(ck_linear,sample) eq 0;
assert Evaluate(ck_cubic,sample) eq 0;
assert Evaluate(boundary,sample) ne 0;
sample01 := EvalList(orbit01,sample);
sample12 := EvalList(orbit12,sample);
rank01, matrix01, primes01 := SquareclassRank(sample01);
rank12, matrix12, primes12 := SquareclassRank(sample12);
assert rank01 eq 4;
assert rank12 eq 4;

// Exact numerical form of the general criterion for an even monic
// split sextic.  The pair (A,B) denotes the class B-A.
function ExactHalvingCriterion(roots, iA, iB)
    A := Q!roots[iA];
    B := Q!roots[iB];
    remaining := [Q!roots[k] : k in [1..6] | k ne iA and k ne iB];
    D := &*[A-c : c in remaining];
    qs := [ -D*(c-B)/(c-A) : c in remaining ];
    return &and[IsSquare(q) : q in qs], qs;
end function;

function ExactJacobianHalvingCriterion(roots, iA, iB)
    f := &*[x-Q!a : a in roots];
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    T2 := J![(x-Q!roots[iA])*(x-Q!roots[iB]),Q!0];
    assert T2 ne J!0 and 2*T2 eq J!0;
    divisible, half := IsDivisibleBy(T2,2);
    if divisible then
        assert 2*half eq T2 and Order(half) eq 4;
    end if;
    return divisible;
end function;

sample_roots := [Q!0] cat [r^2 : r in sample];
for i in [1..5] do
    for j in [i+1..6] do
        criterion, qs := ExactHalvingCriterion(sample_roots,i,j);
        magma_answer := ExactJacobianHalvingCriterion(sample_roots,i,j);
        assert criterion eq magma_answer;
    end for;
end for;

// A positive control catches the omitted-common-squareclass error that
// arises if one tests only cross-ratio quotients.
positive_roots := [Q!-12,Q!-11,Q!-10,Q!-4,Q!-3,Q!-2];
positive, positive_qs := ExactHalvingCriterion(positive_roots,1,6);
assert positive;
assert positive_qs eq [Q!1296,Q!576,Q!36,Q!16];
assert ExactJacobianHalvingCriterion(positive_roots,1,6);

// A negative control where all cross-ratios have the same squareclass,
// but the common Stoll/Zarhin squareclass is nontrivial.
// Here the four cross-ratios are 2,8,18,32, hence all have the
// same squareclass, but their common exact radicand class is not a
// square.  Integral roots keep Magma's genus-2 divisibility routine on
// its certified integral-model path.
negative_roots := [Q!0,Q!3689,Q!-3689,Q!-527,Q!-217,Q!-119];
negative, negative_qs := ExactHalvingCriterion(negative_roots,1,2);
assert not negative;
assert not ExactJacobianHalvingCriterion(negative_roots,1,2);

print "ELKIES22210_HALVING_COVERS";
print "CK_EQUATIONS", ck_linear, ck_cubic;
print "SMOOTH_BOUNDARY_PRODUCT",
      "prod_i(r_i) * prod_{i<j}(r_i^2-r_j^2)";
print "ORBIT_01_RADICANDS", orbit01;
print "ORBIT_12_RADICANDS_FACTORED",
      "[-(a1-a3)(a1-a4)(a1-a5),",
      " (a3-a2)(a1-a4)(a1-a5),",
      " (a4-a2)(a1-a3)(a1-a5),",
      " (a5-a2)(a1-a3)(a1-a4)]";
print "SAMPLE", sample;
print "SAMPLE_01_VALUES", sample01, "SQUARECLASS_RANK", rank01,
      "GENERIC_DEGREE", 2^rank01;
print "SAMPLE_12_VALUES", sample12, "SQUARECLASS_RANK", rank12,
      "GENERIC_DEGREE", 2^rank12;
print "POSITIVE_CONTROL_RADICANDS", positive_qs;
print "NEGATIVE_CONTROL_RADICANDS", negative_qs;
print "ORBIT_01_UNIT_CIRCLE_CHART_VERIFIED", true;
print "ORBIT_01_DIRECT_EQUATIONS",
      "sum_j 1/(1+t_j^2)=3/2, prod_j(1-t_j^2)=16";
print "ORBIT_01_SQUARE_SYMMETRIC_EQUATIONS",
      "e4(t_j^2)=e1(t_j^2)+5, e3(t_j^2)=e2(t_j^2)-10";
print "CK_RATIONAL_TWO_PARAMETER_CHART_AND_INVERSE_VERIFIED", true;
print "ALL_15_ELKIES_CLASSES_MATCH_MAGMA", true;
print "ZARHIN_HALF_DOUBLES_EXACTLY", true;
print "DONE";
quit;
