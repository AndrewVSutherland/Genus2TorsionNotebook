// [2,22] order-11 subsurface geometry probe, STEP 1: fiber degree.
//
// At a base (r,s) the residual resultant Rd(d) = Res_q(N0/G, N1/G) has
// degree ~1529; its roots split among the order-5 subcomponent, the
// order-11 subcomponent S11, degenerate boundary components, etc.
// Over F_p we factor Rd completely and TAG each irreducible factor by the
// CF order of the D_inf class at its points (roots taken in F_{p^k},
// k <= 3; CF and the curve construction are field-generic).  The total
// degree mass tagged 11 is the fiber degree n11 of S11 over the
// (r,s)-plane -- it must be CONSTANT across generic bases and primes.
// This n11 drives step 2 (interpolation of the S11 equation + genus).
//
// Run: magma -b code/blp22_ord11_geom_step1.m > results/blp22_ord11_geom_step1.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(4*10^9);

load "code/blp22_locus_rm_test2.m0";   // DeriveT

function SqrtPolyPart6(f)
    P := Parent(f); xx := P.1; sp := xx^3;
    for k in [1..3] do
        dd := f - sp^2;
        if Degree(dd) le 2 then break; end if;
        sp := sp + (Coefficient(dd, 6-k)/(2*Coefficient(sp, 3)))*xx^(3-k);
    end for;
    return sp;
end function;
function CFOrderB(f, maxsteps, maxord)
    P := Parent(f); sp := SqrtPolyPart6(f);
    Pi := P!0; Qi := P!1; total := 0;
    for i in [0..maxsteps] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + sp) div Qi;
        total +:= Degree(ai);
        if total gt maxord then return 0; end if;
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;

// tag one irreducible factor h(d) over Fp: return CF order at a point above
// a root of h, or -1 (degenerate/undecided)
function TagFactor(h, N0r, N1r, r, sv, p)
    k := Degree(h);
    if k gt 60 then return -2; end if;  // tag everything up to deg 60 (CF-only, cheap over F_{p^k}); residual -2 mass = parasitic giants
    Fq := ext< GF(p) | k >;
    rts := Roots(PolynomialRing(Fq)!h);
    if #rts eq 0 then return -1; end if;
    d0 := rts[1][1];
    if Fq!sv^2 + d0 eq 0 then return -1; end if;
    PE := PolynomialRing(Fq);
    // specialize the bivariate polys at d = d0 via coefficient extraction:
    // N0r = sum_{i,j} c_ij q^i d^j  -> A0(q) = sum_i (sum_j c_ij d0^j) q^i
    dq0 := Degree(N0r, 1); dq1 := Degree(N1r, 1);
    A0 := PE!0; A1 := PE!0;
    for i in [0..dq0] do
        co := Fq!0;
        cf := Coefficient(N0r, 1, i);   // poly in d
        for j in [0..Degree(cf, 2)] do
            co +:= (Fq!Coefficient(cf, 2, j))*d0^j;
        end for;
        A0 +:= co*PE.1^i;
    end for;
    for i in [0..dq1] do
        co := Fq!0;
        cf := Coefficient(N1r, 1, i);
        for j in [0..Degree(cf, 2)] do
            co +:= (Fq!Coefficient(cf, 2, j))*d0^j;
        end for;
        A1 +:= co*PE.1^i;
    end for;
    if A0 eq 0 or A1 eq 0 then return -1; end if;
    gg := GCD(A0, A1);
    if Degree(gg) lt 1 then return -1; end if;
    // take a root of gg over an extension if needed (degree cap 2)
    fac := Factorization(gg);
    q0 := Fq!0; found := false;
    for fp in fac do
        if Degree(fp[1]) eq 1 then
            q0 := -Coefficient(fp[1],0)/Coefficient(fp[1],1); found := true; break;
        end if;
    end for;
    if not found then
        for fp in fac do
            if Degree(fp[1]) eq 2 then
                Fq2 := ext< Fq | 2 >;
                rr := Roots(PolynomialRing(Fq2)!fp[1]);
                if #rr gt 0 then
                    // redo in Fq2
                    q02 := rr[1][1];
                    PE2 := PolynomialRing(Fq2); X2 := PE2.1;
                    rE := Fq2!GF(p)!r; sE := Fq2!GF(p)!sv; dE := Fq2!d0;
                    den0 := 2*(sE^2 + dE) - sE*(sE - rE);
                    if den0 eq 0 then return -1; end if;
                    p0 := (2*(rE-1)*(sE^2+dE) + (sE-rE)*(sE^2+q02))/den0;
                    c0 := (rE-1-p0)/2;
                    if c0 eq 0 then return -1; end if;
                    g0 := (X2-rE)*(X2^2+p0*X2+q02);
                    f := g0*(g0 + 4*c0*(X2^2+dE));
                    if Degree(f) ne 6 or Degree(GCD(f,Derivative(f))) gt 0 then return -1; end if;
                    return CFOrderB(f, 60, 45);
                end if;
            end if;
        end for;
        return -1;
    end if;
    PEX := PolynomialRing(Fq); X := PEX.1;
    rE := Fq!GF(p)!r; sE := Fq!GF(p)!sv;
    den0 := 2*(sE^2 + d0) - sE*(sE - rE);
    if den0 eq 0 then return -1; end if;
    p0 := (2*(rE-1)*(sE^2+d0) + (sE-rE)*(sE^2+q0))/den0;
    c0 := (rE-1-p0)/2;
    if c0 eq 0 then return -1; end if;
    g0 := (X-rE)*(X^2+p0*X+q0);
    f := g0*(g0 + 4*c0*(X^2+d0));
    if Degree(f) ne 6 or Degree(GCD(f,Derivative(f))) gt 0 then return -1; end if;
    return CFOrderB(f, 60, 45);
end function;

for p in [499] do
    Fp := GF(p);
    K2<qv, dv> := RationalFunctionField(Fp, 2);
    P2 := PolynomialRing(Fp, 2);
    for base in [1..2] do
        r := Fp!Random(Fp); sv := Fp!Random(Fp);
        if sv eq r then continue; end if;
        rK := K2!r; sK := K2!sv;
        den := 2*(sK^2 + dv) - sK*(sK - rK);
        if den eq 0 then continue; end if;
        pf := (2*(rK-1)*(sK^2+dv) + (sK-rK)*(sK^2+qv)) / den;
        cf := (rK - 1 - pf)/2;
        t0, t1, ok := DeriveT(qv - rK*pf, -rK*qv + 2*cf*dv, cf, K2!dv);
        if not ok or t0 eq 0 or t1 eq 0 then continue; end if;
        N0 := P2!Numerator(t0); N1 := P2!Numerator(t1);
        G := GCD(N0, N1);
        N0r := N0 div G; N1r := N1 div G;
        Rd := Resultant(N0r, N1r, P2.1);
        if Rd eq 0 then continue; end if;
        Rdn := UnivariatePolynomial(Rd);
        Rdn := PolynomialRing(Fp)!Rdn;
        fac := SquarefreeFactorization(Rdn);
        mass := AssociativeArray();
        for pr in fac do
            for fp in Factorization(pr[1]) do
                tag := TagFactor(fp[1], N0r, N1r, r, sv, p);
                key := tag;
                dm := Degree(fp[1]);
                if IsDefined(mass, key) then mass[key] +:= dm; else mass[key] := dm; end if;
            end for;
        end for;
        printf "BASE p=%o r=%o s=%o deg(Rd)=%o MASS:", p, r, sv, Degree(Rdn);
        for k in Sort(Setseq(Keys(mass))) do printf " %o:%o", k, mass[k]; end for;
        printf "\n";
    end for;
end for;
printf "GEOM_STEP1_DONE\n";
quit;
