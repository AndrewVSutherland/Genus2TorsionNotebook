//////////////////////////////////////////////////////////////////////
//  Constructive [4,12] search from the reduced s=m^2 condition.
//
//  For each rational (r,z), find rational roots u of Q4.  For each
//  resulting independent 2-torsion class beta, specialize the generic
//  halving condition to a univariate polynomial F_s(s).  Rational
//  square roots s=m^2 give explicit m,n,A,B and hence a Mumford half
//      H = [ X^2 + A X + B, (X-beta)(mX+n) mod (X^2 + A X + B) ].
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);
PX<X> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function RationalSquareRoot(q)
    if q lt 0 then
        return false, Q!0;
    end if;
    num := Integers()!Numerator(q);
    den := Integers()!Denominator(q);
    if not IsSquare(num) or not IsSquare(den) then
        return false, Q!0;
    end if;
    _, rn := IsSquare(num);
    _, rd := IsSquare(den);
    return true, Q!rn/rd;
end function;

function PrimitiveUnivariate(f)
    R := Parent(f);
    coeffs := Coefficients(f);
    if #coeffs eq 0 then
        return f;
    end if;
    dens := [Denominator(c) : c in coeffs | c ne 0];
    if #dens gt 0 then
        L := LCM(dens);
        f := R!(L*f);
    end if;
    coeffs := Coefficients(f);
    nums := [Integers()!c : c in coeffs | c ne 0];
    if #nums eq 0 then
        return f;
    end if;
    content := GCD([Abs(c) : c in nums]);
    if content gt 1 then
        f := R!(f/content);
    end if;
    return f;
end function;

function M12Data(a, r)
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    W := h^2 + 4*a*x^2*T*(T+1);
    Q4 := ExactQuotient(W, T+1);
    return W, T, h, Q4;
end function;

function OddQuinticAtRoot(W, w)
    out := PX!0;
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
    return PX!(L^2*f), L;
end function;

function IrreducibleFrobeniusCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67] do
        try
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, PX!0;
end function;

function SConditionData(f5, beta)
    Ps<S> := PolynomialRing(Q);
    PXS<XS> := PolynomialRing(FieldOfFractions(Ps));

    fS := PXS!0;
    for i in [0..Degree(f5)] do
        fS +:= Coefficient(f5, i)*XS^i;
    end for;
    g, rem := Quotrem(fS, XS - Ps!beta);
    if rem ne 0 then
        return false, Ps!0, Ps!0, Ps!0, Ps!0, Ps!0, Ps!0, Ps!0;
    end if;

    g4 := Coefficient(g, 4);
    g3 := Coefficient(g, 3);
    g2 := Coefficient(g, 2);
    g1 := Coefficient(g, 1);
    g0 := Coefficient(g, 0);

    A := (g3 - S)/(2*g4);
    D := g2 + beta*S - g4*A^2;
    Pfun := 2*(beta + A);
    Q0 := g1 - A*D;
    L := 8*(g4*beta - S)*(beta + A) + 4*D;
    Cfun := 4*(g4*beta - S)*Q0 + 4*g4*g0 - D^2;
    F := Cfun^2 + S*Pfun*Cfun*L - S*Q0*L^2;
    return true, PrimitiveUnivariate(Ps!F), Ps!A, Ps!D, Ps!Pfun, Ps!Q0,
           Ps!L, Ps!Cfun;
end function;

function ConstructHalfDataFromT(f5, beta, s0, t0, A0, D0, P0, Q00)
    N0 := P0*t0 + Q00;

    okm, m0 := RationalSquareRoot(s0);
    if not okm then
        return false, Q!0, Q!0, PX!0, PX!0;
    end if;

    if m0 eq 0 then
        if t0 ne 0 then
            return false, Q!0, Q!0, PX!0, PX!0;
        end if;
        okn, n0 := RationalSquareRoot(N0);
        if not okn then
            return false, Q!0, Q!0, PX!0, PX!0;
        end if;
    else
        n0 := t0/m0;
        if n0^2 ne N0 then
            return false, Q!0, Q!0, PX!0, PX!0;
        end if;
    end if;

    g := ExactQuotient(f5, X-beta);
    g4 := Coefficient(g, 4);
    B0 := (D0 - 2*t0)/(2*g4);
    U := X^2 + A0*X + B0;
    Vfull := (X-beta)*(m0*X+n0);
    V := Vfull mod U;
    E := g - (X-beta)*(m0*X+n0)^2 - g4*U^2;
    if E ne 0 then
        return false, Q!0, Q!0, PX!0, PX!0;
    end if;
    return true, m0, n0, U, V;
end function;

function ConstructHalfData(f5, beta, s0, A0, D0, P0, Q00, L0, C0)
    if L0 eq 0 then
        return false, Q!0, Q!0, PX!0, PX!0;
    end if;
    t0 := -C0/L0;
    return ConstructHalfDataFromT(f5, beta, s0, t0, A0, D0, P0, Q00);
end function;

print "M(12) full extra-Weierstrass constructive s=m^2 [4,12] search";
print "height", height;

params := RationalParametersOfHeight(height);
hits := [];
checked_split := 0;
extra_roots := 0;
s_polys := 0;
s_roots := 0;
square_s := 0;
zero_s_polys := 0;
special_s_roots := 0;
special_square_s := 0;
constructed := 0;
verified_halves := 0;
order12_tests := 0;

for r in params do
    for z in params do
        if #hits ge max_hits then
            break r;
        end if;
        if r eq -1 or z^2 eq 1 or z eq 0 then
            continue;
        end if;

        a := (1-z^2)/(4*(r+1));
        if a eq 0 then
            continue;
        end if;

        W, T, h, Q4 := M12Data(a, r);
        if Degree(W) ne 6 or Discriminant(W) eq 0 then
            continue;
        end if;

        rootsT := Roots(T+1);
        if #rootsT lt 2 then
            continue;
        end if;
        checked_split +:= 1;

        rootsQ := [ rt[1] : rt in Roots(Q4) | rt[2] eq 1 ];
        if #rootsQ eq 0 then
            continue;
        end if;
        extra_roots +:= #rootsQ;

        for wd in rootsT do
            w := wd[1];
            if w eq 0 then
                continue;
            end if;

            f5 := OddQuinticAtRoot(W, w);
            if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
                continue;
            end if;

            fI, Lscale := IntegralModelPolynomial(f5);
            if Discriminant(fI) eq 0 then
                continue;
            end if;
            CI := HyperellipticCurve(fI);
            JI := Jacobian(CI);

            Y0 := Evaluate(h, 0);
            Xp := -1/w;
            Yp := Y0*Xp^3;
            DI := JI![X-Xp, Lscale*Yp];
            try
                ordD := Order(DI);
            catch e
                continue;
            end try;
            order12_tests +:= 1;
            if ordD ne 12 then
                continue;
            end if;
            TdivI := 6*DI;

            Ccurve := HyperellipticCurve(f5);
            J := Jacobian(Ccurve);

            for u in rootsQ do
                if u eq w then
                    continue;
                end if;
                beta := 1/(u-w);

                TbetaI := JI![X-beta, Q!0];
                if TbetaI eq JI!0 or TbetaI eq TdivI then
                    continue;
                end if;

                ok, Fpoly, Apoly, Dpoly, Ppoly, Q0poly, Lpoly, Cpoly :=
                    SConditionData(f5, beta);
                if not ok then
                    continue;
                end if;
                s_polys +:= 1;

                if Fpoly eq 0 then
                    zero_s_polys +:= 1;
                else
                    for rt in Roots(Fpoly) do
                        s0 := rt[1];
                    s_roots +:= 1;
                    if Evaluate(Lpoly, s0) eq 0 then
                        continue;
                    end if;
                    okm := RationalSquareRoot(s0);
                    if not okm then
                        continue;
                    end if;
                    square_s +:= 1;

                    A0 := Evaluate(Apoly, s0);
                    D0 := Evaluate(Dpoly, s0);
                    P0 := Evaluate(Ppoly, s0);
                    Q00 := Evaluate(Q0poly, s0);
                    L0 := Evaluate(Lpoly, s0);
                    C0 := Evaluate(Cpoly, s0);

                    okH, m0, n0, U, V := ConstructHalfData(f5, beta, s0, A0, D0, P0, Q00, L0, C0);
                    if not okH then
                        continue;
                    end if;
                    constructed +:= 1;

                    H := J![U, V];
                    Tbeta := J![X-beta, Q!0];
                    if 2*H ne Tbeta and 2*H ne -Tbeta then
                        continue;
                    end if;
                    verified_halves +:= 1;

                    simple, pcert, Lp := IrreducibleFrobeniusCertificate(f5);
                    Append(~hits, <r,z,a,w,u,beta,s0,m0,n0,U,V,f5,simple,pcert,Lp>);
                    print "HIT";
                    print "  r,z,a", r, z, a;
                    print "  w,u,beta", w, u, beta;
                    print "  s,m,n", s0, m0, n0;
                    print "  U,V", U, V;
                    print "  simple", simple, "prime", pcert, "L", Lp;
                    print "  f5", f5;
                    end for;
                end if;

                Gspecial := GCD(PrimitiveUnivariate(Lpoly), PrimitiveUnivariate(Cpoly));
                if Gspecial ne 0 and Degree(Gspecial) ge 1 then
                    for srt in Roots(Gspecial) do
                        s0 := srt[1];
                        special_s_roots +:= 1;
                        okm := RationalSquareRoot(s0);
                        if not okm then
                            continue;
                        end if;
                        special_square_s +:= 1;

                        A0 := Evaluate(Apoly, s0);
                        D0 := Evaluate(Dpoly, s0);
                        P0 := Evaluate(Ppoly, s0);
                        Q00 := Evaluate(Q0poly, s0);
                        Pt<Tt> := PolynomialRing(Q);
                        tpoly := Tt^2 - s0*P0*Tt - s0*Q00;
                        for trt in Roots(tpoly) do
                            t0 := trt[1];
                            okH, m0, n0, U, V := ConstructHalfDataFromT(f5, beta, s0, t0, A0, D0, P0, Q00);
                            if not okH then
                                continue;
                            end if;
                            constructed +:= 1;

                            H := J![U, V];
                            Tbeta := J![X-beta, Q!0];
                            if 2*H ne Tbeta and 2*H ne -Tbeta then
                                continue;
                            end if;
                            verified_halves +:= 1;

                            simple, pcert, Lp := IrreducibleFrobeniusCertificate(f5);
                            Append(~hits, <r,z,a,w,u,beta,s0,m0,n0,U,V,f5,simple,pcert,Lp>);
                            print "HIT special";
                            print "  r,z,a", r, z, a;
                            print "  w,u,beta", w, u, beta;
                            print "  s,m,n,t", s0, m0, n0, t0;
                            print "  U,V", U, V;
                            print "  simple", simple, "prime", pcert, "L", Lp;
                            print "  f5", f5;
                        end for;
                    end for;
                end if;
            end for;
        end for;

        if progress_interval gt 0 and checked_split mod progress_interval eq 0 then
            print "checked_split", checked_split, "extra_roots", extra_roots,
                  "s_polys", s_polys, "s_roots", s_roots,
                  "square_s", square_s, "zero_s_polys", zero_s_polys,
                  "special_s_roots", special_s_roots,
                  "special_square_s", special_square_s,
                  "constructed", constructed,
                  "verified", verified_halves, "hits", #hits;
        end if;
    end for;
end for;

print "Done";
print "checked_split", checked_split;
print "extra_roots", extra_roots;
print "order12_tests", order12_tests;
print "s_polys", s_polys;
print "s_roots", s_roots;
print "square_s", square_s;
print "zero_s_polys", zero_s_polys;
print "special_s_roots", special_s_roots;
print "special_square_s", special_square_s;
print "constructed", constructed;
print "verified_halves", verified_halves;
print "hits", #hits;
for H in hits do
    print H;
end for;

quit;
