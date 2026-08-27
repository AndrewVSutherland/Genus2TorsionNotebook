// DB-driven reducible-2-torsion gluing sweep (HLP Prop 4 / BHLS A:glue2).
// Pool: LMFDB curves with torsion [8],[10],[12],[2,4],[2,6],[2,8].
// Buckets: full-split 2-torsion (always compatible) and partial (match the
// quadratic 2-division field).  Prefilter pairs by
//   gcd_p ( #E(F_p) * #F(F_p) ) >= MinOrd   (cached elliptic counts),
// then glue (all matchings) and compute exact torsion.
// Params: Sh, NSh, MinOrd.
SetColumns(0);
SetSeed(1);
if not assigned Sh then Sh := 0; elif Type(Sh) eq MonStgElt then Sh := StringToInteger(Sh); end if;
if not assigned NSh then NSh := 3; elif Type(NSh) eq MonStgElt then NSh := StringToInteger(NSh); end if;
if not assigned MinOrd then MinOrd := 48; elif Type(MinOrd) eq MonStgElt then MinOrd := StringToInteger(MinOrd); end if;
SetMemoryLimit(6*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

// load pool via python-written json -> convert offline: read simple lines file instead
pool := [];  // <label, ainvs, torsstring>
lines := Split(Read("data/split_ec_torsion_pool.txt"), "\n");
for L in lines do
    if #L eq 0 then continue; end if;
    parts := Split(L, "|");
    ai := eval parts[2];
    Append(~pool, <parts[1], ai, parts[3]>);
end for;
printf "POOL %o curves\n", #pool;

ps := [ p : p in [11,13,17,19,23,29,31,37,41,43] ];
data := [];
for rec in pool do
    E := EllipticCurve([Q!a : a in rec[2]]);
    W := WeierstrassModel(E);
    ai := aInvariants(W);
    f := x^3 + ai[4]*x + ai[5];
    fac := Factorization(f);
    degs := Sort([Degree(t[1]) : t in fac]);
    if degs eq [1,1,1] then
        bucket := 0;  // full split
    elif degs eq [1,2] then
        q2 := [t[1] : t in fac | Degree(t[1]) eq 2][1];
        d2 := Discriminant(q2);
        bucket := SquarefreeFactorization(Numerator(d2)*Denominator(d2));
        if bucket eq 1 then bucket := 0; end if;   // splits after all
    else
        continue;  // irreducible: handled in earlier campaign
    end if;
    counts := [];
    okp := true;
    for p in ps do
        if Numerator(Discriminant(E)) mod p eq 0 then Append(~counts, 0);
        else Append(~counts, #ChangeRing(E, GF(p))); end if;
    end for;
    Append(~data, <rec[1], f, bucket, counts, rec[3]>);
end for;
printf "PREPARED %o (buckets computed)\n", #data;

function GlueP4(f, g, K)
    rf := [ r[1] : r in Roots(PolynomialRing(K)!f) ];
    rg0 := [ r[1] : r in Roots(PolynomialRing(K)!g) ];
    if #rf ne 3 or #rg0 ne 3 then return []; end if;
    Df := Discriminant(f); Dg := Discriminant(g);
    out := [];
    for s in SymmetricGroup(3) do
        rg := [ rg0[1^s], rg0[2^s], rg0[3^s] ];
        al := rf; be := rg;
        d_al := [ al[3]-al[2], al[1]-al[3], al[2]-al[1] ];
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
        ok := true; cq := [];
        for co in cf do
            fl, cr := IsCoercible(Q, co);
            if not fl then ok := false; break; end if;
            Append(~cq, cr);
        end for;
        if not ok then continue; end if;
        hQ := P!cq;
        if Degree(hQ) eq 6 and IsSquarefree(hQ) then Append(~out, hQ); end if;
    end for;
    return out;
end function;

seen := {};
nesc := 0;
for i in [1..#data] do
    if (i mod NSh) ne Sh then continue; end if;
    di := data[i];
    for j in [i+1..#data] do
        dj := data[j];
        if di[3] ne dj[3] then continue; end if;   // bucket match
        // prefilter: gcd of products
        g := 0;
        for k in [1..#ps] do
            if di[4][k] eq 0 or dj[4][k] eq 0 then continue; end if;
            g := GCD(g, di[4][k]*dj[4][k]);
            if g ne 0 and g lt MinOrd then break; end if;
        end for;
        if g eq 0 or g lt MinOrd then continue; end if;
        f := di[2]; gg := dj[2];
        K := di[3] eq 0 select Q else NumberField(Polynomial(Q,[-di[3],0,1]));
        if Type(K) eq FldRat then
            hs := GlueP4(f, gg, Q);
        else
            hs := GlueP4(f, gg, K);
        end if;
        for hQ in hs do
            ih := Hash(G2Invariants(HyperellipticCurve(hQ)));
            if ih in seen then continue; end if;
            Include(~seen, ih);
            nesc +:= 1;
            C := HyperellipticCurve(hQ);
            Cm := C;
            try Cm := ReducedMinimalWeierstrassModel(C); catch e ; end try;
            T := TorsionSubgroup(Jacobian(SimplifiedModel(Cm)));
            printf "DBGLUE %o(%o) x %o(%o) gcdbound %o TORSION %o order %o\n",
                di[1], di[5], dj[1], dj[5], g, Invariants(T), #T;
        end for;
    end for;
    if i mod 50 eq 0 then printf "PROGRESS i=%o/%o escalated %o\n", i, #data, nesc; end if;
end for;
printf "DBGLUE_DONE shard %o/%o\n", Sh, NSh;
quit;
