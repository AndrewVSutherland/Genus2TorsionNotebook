// Systematic HLP-Prop-4 gluing sweep: pairs of elliptic curves with
// N-torsion (N in {5,7,9}), trivial rational 2-torsion, ISOMORPHIC
// 2-division fields, glued along 2-torsion via the explicit model of
// HLP Prop 4.  Targets: [5,5], [9,9] (in nobody's table), second [63],
// more [7,7]/[35]/[45].
// Validation: must rediscover HLP (7,9;-16/3,4)->63 and (7,7;7,-14/13)->[7,7].
// Params: Types (string like "55" / "77_57" / "99_59_79"), Hnum, Hden.
SetColumns(0);
SetSeed(1);
if not assigned Types then Types := "77_79"; end if;
if not assigned Hnum then Hnum := 16; elif Type(Hnum) eq MonStgElt then Hnum := StringToInteger(Hnum); end if;
if not assigned Hden then Hden := 13; elif Type(Hden) eq MonStgElt then Hden := StringToInteger(Hden); end if;
SetMemoryLimit(4*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

// Tate normal form from Kubert (b,c): y^2+(1-c)xy-by = x^3-bx^2
function TateE(b, c)
    return EllipticCurve([1-c, -b, -b, 0, 0]);
end function;
function KubertBC(N, d)
    // returns b, c, ok
    if N eq 5 then return d, d, true; end if;
    if N eq 7 then return d^3-d^2, d^2-d, true; end if;
    if N eq 9 then
        c := d^2*(d-1);
        b := c*(d*(d-1)+1);
        return b, c, true;
    end if;
    return 0, 0, false;
end function;
// validate the parametrizations on known d
for N in [5,7,9] do
    okN := false;
    for d0 in [Q!2, Q!3, Q!(1/2)] do
        b0, c0, _ := KubertBC(N, d0);
        try
            E0 := TateE(b0, c0);
            if Order(E0![0,0]) eq N then okN := true; break; end if;
        catch e ; end try;
    end for;
    error if not okN, Sprintf("Kubert parametrization for N=%o FAILED validation", N);
    printf "KUBERT_OK N=%o\n", N;
end for;

function MonicCubic(E)
    // short Weierstrass y^2 = x^3 + a4 x + a6
    W := WeierstrassModel(E);
    a := aInvariants(W);
    return x^3 + a[4]*x + a[5];
end function;

// rationals of bounded height
hts := [];
for den in [1..Hden] do
    for num in [-Hnum..Hnum] do
        if num eq 0 then continue; end if;
        if GCD(Abs(num), den) ne 1 then continue; end if;
        Append(~hts, Q!num/den);
    end for;
end for;
printf "PARAMS %o values, types %o\n", #hts, Types;

function GlueP4(f, g, K)
    // all Prop-4 gluings of y^2=f, y^2=g over splitting field K; returns list of sextics in Q[x]
    rf := [ r[1] : r in Roots(PolynomialRing(K)!f) ];
    rg0 := [ r[1] : r in Roots(PolynomialRing(K)!g) ];
    if #rf ne 3 or #rg0 ne 3 then return []; end if;
    Df := Discriminant(f); Dg := Discriminant(g);
    out := [];
    S3 := SymmetricGroup(3);
    for s in S3 do
        rg := [ rg0[1^s], rg0[2^s], rg0[3^s] ];
        al := rf; be := rg;
        d_al := [ al[3]-al[2], al[1]-al[3], al[2]-al[1] ];  // indexed cyclic
        d_be := [ be[3]-be[2], be[1]-be[3], be[2]-be[1] ];
        if 0 in d_al or 0 in d_be then continue; end if;
        a1 := d_al[1]^2/d_be[1] + d_al[3]^2/d_be[3] + d_al[2]^2/d_be[2];
        b1 := d_be[1]^2/d_al[1] + d_be[3]^2/d_al[3] + d_be[2]^2/d_al[2];
        a2 := al[1]*d_be[1] + al[2]*d_be[2] + al[3]*d_be[3];
        b2 := be[1]*d_al[1] + be[2]*d_al[2] + be[3]*d_al[3];
        if a1 eq 0 or a2 eq 0 or b1 eq 0 or b2 eq 0 then continue; end if;
        A := Dg*a1/a2; B := Df*b1/b2;
        PK<xx> := PolynomialRing(K);
        h := -( A*d_al[3]*d_al[2]*xx^2 + B*d_be[3]*d_be[2] )
            *( A*d_al[1]*d_al[3]*xx^2 + B*d_be[1]*d_be[3] )
            *( A*d_al[2]*d_al[1]*xx^2 + B*d_be[2]*d_be[1] );
        cf := Coefficients(h);
        ok := true;
        cq := [];
        for co in cf do
            fl, cr := IsCoercible(Q, co);
            if not fl then ok := false; break; end if;
            Append(~cq, cr);
        end for;
        if not ok then continue; end if;
        hQ := P!cq;
        if Degree(hQ) ne 6 or not IsSquarefree(hQ) then continue; end if;
        Append(~out, hQ);
    end for;
    return out;
end function;

typelist := Split(Types, "_");
seenh := {};
for ty in typelist do
    M := StringToInteger(ty[1]);
    N := StringToInteger(ty[2]);
    printf "TYPE %o x %o\n", M, N;
    for i in [1..#hts] do
        tv := hts[i];
        bM, cM, _ := KubertBC(M, tv);
        okE := true;
        try
            E := TateE(bM, cM);
        catch e okE := false; end try;
        if not okE then continue; end if;
        if Order(E![0,0]) ne M then continue; end if;
        f := MonicCubic(E);
        if not IsIrreducible(f) then continue; end if;
        Df := Discriminant(f);
        jE := jInvariant(E);
        jstart := (M eq N) select i+1 else 1;
        for jdx in [jstart..#hts] do
            uv := hts[jdx];
            bN, cN, _ := KubertBC(N, uv);
            okF := true;
            try
                F := TateE(bN, cN);
            catch e okF := false; end try;
            if not okF then continue; end if;
            if Order(F![0,0]) ne N then continue; end if;
            if jInvariant(F) eq jE then continue; end if;
            g := MonicCubic(F);
            if not IsIrreducible(g) then continue; end if;
            if not IsSquare(Df*Discriminant(g)) then continue; end if;
            Kf := NumberField(f);
            Kg := NumberField(g);
            if not IsIsomorphic(Kf, Kg) then continue; end if;
            printf "PAIRHIT %ox%o t=%o u=%o\n", M, N, tv, uv;
            K := SplittingField(f);
            hs := GlueP4(f, g, K);
            for hQ in hs do
                C := HyperellipticCurve(hQ);
                ii := G2Invariants(C);
                if ii in seenh then continue; end if;
                Include(~seenh, ii);
                Cm := C;
                try Cm := ReducedMinimalWeierstrassModel(C); catch e ; end try;
                T := TorsionSubgroup(Jacobian(SimplifiedModel(Cm)));
                printf "GLUE %ox%o t=%o u=%o TORSION %o order %o\n", M, N, tv, uv, Invariants(T), #T;
            end for;
        end for;
        if i mod 40 eq 0 then printf "PROGRESS type %ox%o i=%o/%o\n", M, N, i, #hts; end if;
    end for;
end for;
printf "GLUE_SWEEP_DONE %o\n", Types;
quit;
