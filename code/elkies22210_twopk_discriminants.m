//////////////////////////////////////////////////////////////////////
//  Doubled-theta loci on the Elkies Clebsch--Klein [2,2,2,10]
//  surface.
//
//  Move the distinguished Weierstrass point x=0 of
//
//      y^2 = x * Product_i (x-r_i^2)
//
//  to infinity.  The resulting odd model is
//
//      Y^2 = Product_i (1-r_i^2 X),       P=(0,1),
//
//  and T=P-infinity is the marked order-5 class.  Up to sign and the
//  S_5 action, the order-10 classes which can possibly be 2R-K are
//
//      A: T  + S_ij,
//      B: 2T + S_0i,
//      C: 2T + S_ij.
//
//  In each case a polynomial ell cuts out the prescribed part of the
//  divisor plus a residual effective divisor of degree 2.  RawU below
//  cuts out that residual divisor.  Hence Disc_X(RawU)=0 is precisely
//  the affine doubled-theta condition.  The residual class is the
//  negative of the displayed order-10 class, which has the same
//  repeated-support condition.
//
//  Typical runs (always use a shell memory cap for factor mode):
//
//      magma -b mode:=derive code/elkies22210_twopk_discriminants.m
//      magma -b mode:=factor case_name:=A \
//          code/elkies22210_twopk_discriminants.m
//      magma -b mode:=verify code/elkies22210_twopk_discriminants.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "derive";
end if;
if not assigned case_name then
    case_name := "all";
end if;

Q := Rationals();
R<t,m> := PolynomialRing(Q, 2);
RX<X> := PolynomialRing(R);

// Full rational chart of the labelled smooth Clebsch--Klein surface.
rr := [
    1+t*(t+2)*m,
    t*m*(m-t-2),
    -1+m+t*(t+1)*m^2,
    1+t-m-t*m^2,
    -(1+t)*(1+t*m^2)
];
assert &+rr eq 0;
assert &+[u^3 : u in rr] eq 0;

aa := [u^2 : u in rr];
L := [1-aa[i]*X : i in [1..5]];
F := &*L;
fprime0 := Coefficient(F, 1);
tangent_slope := fprime0/2;

function ProductExcept(seq, omitted)
    ans := Universe(seq)!1;
    for i in [1..#seq] do
        if i notin omitted then
            ans *:= seq[i];
        end if;
    end for;
    return ans;
end function;

function ExactQuotient(num, den)
    ok, quo := IsDivisibleBy(num, den);
    assert ok;
    return quo;
end function;

// Case A: T+S_ij.  The quadratic ell=L_i*L_j passes through P,W_i,W_j.
// Since ell^2-F has the additional root X=0, the two residual roots are
// cut out by (Product_{k notin {i,j}} L_k-L_i*L_j)/X.
function CaseARawU(i, j)
    pair := L[i]*L[j];
    other := ProductExcept(L, {i,j});
    U := ExactQuotient(other-pair, X);
    assert Degree(U) eq 2;
    return U, pair;
end function;

// Case B: 2T+S_0i.  The quadratic ell is tangent to the curve at P and
// vanishes at W_i.  The remaining two intersections give RawU.
function CaseBRawU(i)
    ell := L[i]*(1+(tangent_slope+aa[i])*X);
    U := ExactQuotient(ell^2-F, X^2*L[i]);
    assert Degree(U) eq 2;
    return U, ell;
end function;

// Case C: 2T+S_ij.  A cubic ell is tangent at P and vanishes at W_i,W_j.
// It has pole order 6, so the four prescribed intersections leave two.
function CaseCRawU(i, j)
    ell := L[i]*L[j]
        *(1+(tangent_slope+aa[i]+aa[j])*X);
    U := ExactQuotient(ell^2-F, X^2*L[i]*L[j]);
    assert Degree(U) eq 2;
    return U, ell;
end function;

function QuadraticDiscriminant(U)
    assert Degree(U) eq 2;
    return Coefficient(U,1)^2
        - 4*Coefficient(U,0)*Coefficient(U,2);
end function;

UA, ellA := CaseARawU(1,2);
UB, ellB := CaseBRawU(1);
UC, ellC := CaseCRawU(1,2);
dA := QuadraticDiscriminant(UA);
dB := QuadraticDiscriminant(UB);
dC := QuadraticDiscriminant(UC);

function Bidegree(g)
    mons := Monomials(g);
    if #mons eq 0 then
        return <-1,-1>;
    end if;
    return <Maximum([Degree(mon,1) : mon in mons]),
            Maximum([Degree(mon,2) : mon in mons])>;
end function;

procedure PrintBasic(label, U, disc)
    print "CASE", label;
    print "RAW_U", U;
    print "DISC_BIDEGREE", Bidegree(disc);
    print "DISC_TOTAL_DEGREE", TotalDegree(disc);
    print "DISC_TERMS", #Terms(disc);
    print "DISC", disc;
end procedure;

procedure FactorOne(label, disc)
    print "FACTOR_START", label,
          "bidegree", Bidegree(disc),
          "total_degree", TotalDegree(disc),
          "terms", #Terms(disc);
    time fac := Factorization(disc);
    print "FACTOR", label, fac;
    print "FACTOR_END", label;
end procedure;

function EvalRX(poly, vals, xvar)
    Qx := Parent(xvar);
    ans := Qx!0;
    for k in [0..Degree(poly)] do
        ans +:= Q!Evaluate(Coefficient(poly,k), vals)*xvar^k;
    end for;
    return ans;
end function;

procedure VerifyAt(vals)
    Qx<x> := PolynomialRing(Q);
    rq := [Q!Evaluate(r,vals) : r in rr];
    aq := [r^2 : r in rq];
    fq := EvalRX(F, vals, x);
    assert Degree(fq) eq 5 and Discriminant(fq) ne 0;

    C := HyperellipticCurve(fq);
    J := Jacobian(C);
    T5 := J![x,Q!1];
    assert Order(T5) eq 5;

    S01 := J![x-1/aq[1],Q!0];
    S12 := J![(x-1/aq[1])*(x-1/aq[2]),Q!0];
    assert S01 ne J!0 and 2*S01 eq J!0;
    assert S12 ne J!0 and 2*S12 eq J!0;

    Udata := [UA,UB,UC];
    Edata := [ellA,ellB,ellC];
    Ddata := [T5+S12,2*T5+S01,2*T5+S12];
    labels := ["A","B","C"];
    for k in [1..3] do
        Uqraw := EvalRX(Udata[k],vals,x);
        Uq := Uqraw/LeadingCoefficient(Uqraw);
        ellq := EvalRX(Edata[k],vals,x);
        _, vq := Quotrem(ellq,Uq);
        residual := J![Uq,vq];
        assert residual eq -Ddata[k];
        assert Order(Ddata[k]) eq 10;
        assert Discriminant(Uqraw) eq
            Q!Evaluate([dA,dB,dC][k],vals);
        print "VERIFY_CASE", labels[k],
              "order", Order(Ddata[k]),
              "raw_u", Uqraw,
              "disc", Discriminant(Uqraw);
    end for;

    print "VERIFY_OK", vals, "r", rq;
end procedure;

if mode eq "derive" then
    if case_name in {"all","A"} then PrintBasic("A",UA,dA); end if;
    if case_name in {"all","B"} then PrintBasic("B",UB,dB); end if;
    if case_name in {"all","C"} then PrintBasic("C",UC,dC); end if;
elif mode eq "factor" then
    if case_name in {"all","A"} then FactorOne("A",dA); end if;
    if case_name in {"all","B"} then FactorOne("B",dB); end if;
    if case_name in {"all","C"} then FactorOne("C",dC); end if;
elif mode eq "verify" then
    // This maps projectively to Elkies's [1,-8,-7,5,9] example.
    VerifyAt([Q!2/7,Q!-14]);
else
    error "mode must be derive, factor, or verify";
end if;
