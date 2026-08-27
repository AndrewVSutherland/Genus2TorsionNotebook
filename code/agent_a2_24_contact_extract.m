
//////////////////////////////////////////////////////////////////////
//  Cubic-contact 3-torsion machinery + validation.
//
//  Claim: a rational order-3 class D=(u,v) (Mumford, u monic deg 2) on
//  J(y^2=f) gives a cubic-contact identity
//        f = h3^2 - kappa * u^3,    h3 = v + u*(s*x + w),  kappa in Q,
//  i.e. u^3 | (f - h3^2) with quotient a constant.  Then q3 := u is the
//  3-torsion "contact conic" and h3 the contact cubic.
//
//  Derivation of (s,w,kappa): u | (f - v^2) since (root of u, v) lies on
//  the curve, so f - v^2 = u*g4 (deg g4 = 4).  Impose f - h3^2 = kappa u^3:
//     f - (v + u(sx+w))^2 = kappa u^3
//     u*[ g4 - 2v(sx+w) - u(sx+w)^2 ] = kappa u^3
//     g4 - 2v(sx+w) - u(sx+w)^2 = kappa u^2   (deg-4 identity: 5 eqns,
//                                              unknowns s,w,kappa).
//  Solve the linear-in-(s,w,kappa)-after-expansion system (it is
//  quadratic in s,w; solve exactly via the 5 coeff eqns, consistent
//  because D really is 3-torsion).
//
//  This script regenerates the two known simple Z/24 curves from their
//  (r,p,t) and validates the identity, printing (h3,q3,kappa).
//
//  Usage: magma -b agent_a2_24_contact_extract.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function A8f(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv; d := e + 2*pv - rv^2; lambda := rv/tv;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d; q := a*x^2 + b*x + c;
    f := q*(Qpoly^2 + q);
    return f;
end function;

function IntModel(fv)
    L := 1;
    for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
    return P!(L^2*fv);
end function;

// given f (deg 6) and a 3-torsion divisor D=(u,v), u monic deg 2, v deg<=1,
// return ok, h3, q3(=u), kappa  with  f = h3^2 - kappa*u^3.
function ExtractContact(f, u, v)
    if Degree(u) ne 2 then return false, P!0, P!0, Q!0; end if;
    rem := (f - v^2) mod u;
    if rem ne 0 then return false, P!0, P!0, Q!0; end if;
    g4 := (f - v^2) div u;               // deg 4
    // unknowns s,w,kappa; identity  g4 - 2 v (s x + w) - u (s x + w)^2 - kappa u^2 = 0
    // treat as poly identity; expand with symbolic s,w,kappa via a poly ring
    R<S,W,K> := PolynomialRing(Q, 3);
    PR := PolynomialRing(R); X := PR.1;
    uR := PR![R!Coefficient(u,i) : i in [0..Degree(u)]];
    vR := PR![R!Coefficient(v,i) : i in [0..Degree(v)]];
    g4R := PR![R!Coefficient(g4,i) : i in [0..Degree(g4)]];
    lin := S*X + W;
    expr := g4R - 2*vR*lin - uR*lin^2 - K*uR^2;
    eqs := [Coefficient(expr, i) : i in [0..4]];
    I := ideal<R | eqs>;
    V := Variety(I);
    if #V eq 0 then return false, P!0, P!0, Q!0; end if;
    sol := V[1];
    sv := Q!sol[1]; wv := Q!sol[2]; kv := Q!sol[3];
    h3 := v + u*(sv*x + wv);
    return true, h3, u, kv;
end function;

// find a rational order-3 divisor on J via the torsion group
function ThreeTorsionDivisor(J)
    T, mp := TorsionSubgroup(J);
    for g in T do
        if Order(g) eq 3 then
            D := mp(g);
            return true, D[1], D[2];   // u, v
        end if;
    end for;
    return false, _, _;
end function;

curves := [ <5, -5/2, -9/2, "A (chi_17)">, <1/3, -1/9, -1, "B (chi_13)"> ];
for cc in curves do
    rv := cc[1]; pv := cc[2]; tv := cc[3]; nm := cc[4];
    f := IntModel(A8f(rv, pv, tv));
    printf "=== curve %o : r=%o p=%o t=%o ===\n", nm, rv, pv, tv;
    J := Jacobian(HyperellipticCurve(f));
    inv := Invariants(TorsionSubgroup(J));
    printf "  torsion = %o\n", inv;
    ok3, u, v := ThreeTorsionDivisor(J);
    if not ok3 then printf "  no rational order-3 divisor found\n"; continue; end if;
    printf "  3-torsion Mumford: u=%o  v=%o\n", u, v;
    ok, h3, q3, kappa := ExtractContact(f, u, v);
    if not ok then printf "  ExtractContact FAILED\n"; continue; end if;
    printf "  q3   = %o\n", q3;
    printf "  h3   = %o\n", h3;
    printf "  kappa= %o\n", kappa;
    check := f - (h3^2 + kappa*q3^3);
    printf "  IDENTITY f - (h3^2 + kappa*q3^3) = %o  %o\n",
        check, check eq 0 select "[VALID]" else "[FAIL]";
    printf "  lead(f)=%o  lead(h3)^2 + kappa = %o\n",
        LeadingCoefficient(f), LeadingCoefficient(h3)^2 + kappa;
end for;
print "DONE";
quit;
