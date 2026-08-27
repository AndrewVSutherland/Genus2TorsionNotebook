//////////////////////////////////////////////////////////////////////
//  General non-d=0 cubic-contact reconnaissance for [2,24].
//
//  This is not a blind rational search.  It checks the open contact
//  door left by the d=0 rank-0 result:
//
//      f = h3^2 + kappa*(x^2 + U*x + V)^3
//
//  on the A(8) chart, with particular attention to the W-split
//  (2-rank >= 2) locus.  The survey is finite-field/intrinsic:
//  for each A(8) curve over F_l, compute the rational 3-torsion in
//  the Jacobian, extract the Mumford u-polynomial q3, and classify
//  q3 as split/irreducible.  This avoids substituting 3 | #J as a
//  proxy when we actually want contact classes.
//
//  Usage:
//      magma -b Primes:="7,11,13" \
//          code/agent_a2_24_general_contact_recon.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned MemGB then
    MemGB := 4;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
SetMemoryLimit(MemGB*10^9);

if not assigned Primes then
    Primes := "7,11,13";
elif Type(Primes) ne MonStgElt then
    Primes := Sprint(Primes);
end if;
PrimeList := [StringToInteger(s) : s in Split(Primes, ",")];

if not assigned MaxCurvesPerPrime then
    MaxCurvesPerPrime := 0;
elif Type(MaxCurvesPerPrime) eq MonStgElt then
    MaxCurvesPerPrime := StringToInteger(MaxCurvesPerPrime);
end if;

if not assigned SampleLimit then
    SampleLimit := 8;
elif Type(SampleLimit) eq MonStgElt then
    SampleLimit := StringToInteger(SampleLimit);
end if;

Z := Integers();
Q := Rationals();
PQ<xQ> := PolynomialRing(Q);

function A8fQ(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv;
    d := e + 2*pv - rv^2;
    lambda := rv/tv;
    u := pv + rv*tv - 2*rv;
    v := e + rv^2 - rv*pv - rv^2*tv + 3*pv*tv - rv*tv^2;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e)
         + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := xQ^2 + d;
    q := a*xQ^2 + b*xQ + c;
    g8 := xQ^2 + u*xQ + v;
    f := q*(Qpoly^2 + q);
    L := rv*xQ + (pv - rv^2);
    ellBase := -(q + Qpoly*L);
    return f, g8, ellBase, d;
end function;

function IntModel(fv)
    L := 1;
    for i in [0..Degree(fv)] do
        L := LCM(L, Denominator(Coefficient(fv, i)));
    end for;
    return PQ!(L^2*fv), L;
end function;

function ExtractContactQ(f, u, v)
    if Degree(u) ne 2 then
        return false, PQ!0, PQ!0, Q!0;
    end if;
    if (f - v^2) mod u ne 0 then
        return false, PQ!0, PQ!0, Q!0;
    end if;
    g4 := (f - v^2) div u;
    R<S,W,K> := PolynomialRing(Q, 3);
    PR<X> := PolynomialRing(R);
    uR := PR![R!Coefficient(u,i) : i in [0..Degree(u)]];
    vR := PR![R!Coefficient(v,i) : i in [0..Degree(v)]];
    g4R := PR![R!Coefficient(g4,i) : i in [0..Degree(g4)]];
    lin := S*X + W;
    expr := g4R - 2*vR*lin - uR*lin^2 - K*uR^2;
    eqs := [Coefficient(expr, i) : i in [0..4]];
    V := Variety(ideal<R | eqs>);
    if #V eq 0 then
        return false, PQ!0, PQ!0, Q!0;
    end if;
    sol := V[1];
    sv := Q!sol[1];
    wv := Q!sol[2];
    kv := Q!sol[3];
    h3 := v + u*(sv*xQ + wv);
    return true, h3, u, kv;
end function;

function ThreeTorsionDivisorQ(J)
    T, mp := TorsionSubgroup(J);
    for g in T do
        if Order(g) eq 3 then
            D := mp(g);
            return true, D[1], D[2];
        end if;
    end for;
    return false, PQ!0, PQ!0;
end function;

function A8fFF(x, rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv;
    d := e + 2*pv - rv^2;
    lambda := rv/tv;
    u := pv + rv*tv - 2*rv;
    v := e + rv^2 - rv*pv - rv^2*tv + 3*pv*tv - rv*tv^2;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e)
         + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d;
    q := a*x^2 + b*x + c;
    g8 := x^2 + u*x + v;
    f := q*(Qpoly^2 + q);
    L := rv*x + (pv - rv^2);
    ellBase := -(q + Qpoly*L);
    return f, g8, ellBase, d;
end function;

function TwoRankFF(fp)
    degs := [Degree(g[1]) : g in Factorization(fp)];
    k := #degs;
    even := 0;
    for mask in [0..2^k-1] do
        ss := 0;
        for i in [1..k] do
            if (mask div 2^(i-1)) mod 2 eq 1 then
                ss +:= degs[i];
            end if;
        end for;
        if ss mod 2 eq 0 then
            even +:= 1;
        end if;
    end for;
    return Ilog2(even) - 1;
end function;

function IsEvenSextic(fp)
    return Coefficient(fp,1) eq 0 and Coefficient(fp,3) eq 0
           and Coefficient(fp,5) eq 0;
end function;

function ContactClassesFF(J, f)
    // Return one record per nonzero rational order-3 element with a
    // squarefree degree-2 Mumford u-polynomial.  Magma does not implement
    // TorsionSubgroup for finite-field genus-2 Jacobians, so enumerate the
    // Mumford representatives directly.
    P := Parent(f);
    x := P.1;
    F := BaseRing(P);
    O := J!0;
    rows := [];
    seen := AssociativeArray();
    for U in F do
    for V in F do
        u := x^2 + U*x + V;
        disc := Discriminant(u);
        if disc eq 0 then
            continue;
        end if;
        for A in F do
        for B in F do
            v := A*x + B;
            if (f - v^2) mod u ne 0 then
                continue;
            end if;
            ok := true;
            try
                D := elt<J | u, v>;
            catch e
                ok := false;
            end try;
            if not ok or D eq O or 3*D ne O then
                continue;
            end if;
            // D and -D have the same contact conic q3=u; count the conic
            // once, since split/irreducible is a property of q3.
            key := Sprint(u);
            if IsDefined(seen, key) then
                continue;
            end if;
            seen[key] := true;
            split := IsSquare(disc);
            Append(~rows, <u, split>);
        end for;
        end for;
    end for;
    end for;
    return rows;
end function;

procedure Inc(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

procedure AddSample(~samples, item, lim)
    if #samples lt lim then
        Append(~samples, item);
    end if;
end procedure;

print "GENERAL_NON_D0_CONTACT_RECON_START";
print "Primes", PrimeList, "MaxCurvesPerPrime", MaxCurvesPerPrime,
      "SampleLimit", SampleLimit;

print "";
print "Q_ANCHOR_VALIDATION";
Anchors := [
    <Q!5, Q!(-5)/2, Q!(-9)/2, "A_non_d0_irreducible_q3">,
    <Q!1/3, Q!(-1)/9, Q!(-1), "B_d0_split_q3">
];
for A in Anchors do
    rv := A[1]; pv := A[2]; tv := A[3]; name := A[4];
    f, g8, ellBase, dval := A8fQ(rv, pv, tv);
    fInt, Lden := IntModel(f);
    J := Jacobian(HyperellipticCurve(fInt));
    inv := Invariants(TorsionSubgroup(J));
    ok3, u, v := ThreeTorsionDivisorQ(J);
    printf " anchor %o r=%o p=%o t=%o d=%o torsion=%o order3=%o\n",
        name, rv, pv, tv, dval, inv, ok3;
    if ok3 then
        okc, h3, q3, kap := ExtractContactQ(fInt, u, v);
        disc := Discriminant(q3);
        printf "  q3=%o disc=%o split_over_Q=%o contact_identity=%o kappa=%o\n",
            q3, disc, IsSquare(disc), okc, kap;
    end if;
end for;

print "";
print "FINITE_FIELD_W_SPLIT_CONTACT_SURVEY";
print "columns: ell base smooth order8 wsplit wsplit_j3 contact_bases contact_classes split_classes irreducible_classes non_d0_contact_bases non_d0_irred_bases d0_contact_bases";

for ell in PrimeList do
    F := GF(ell);
    P<x> := PolynomialRing(F);
    base := 0;
    smooth := 0;
    order8 := 0;
    wsplit := 0;
    wsplit_j3 := 0;
    contact_bases := 0;
    contact_classes := 0;
    split_classes := 0;
    irreducible_classes := 0;
    non_d0_contact_bases := 0;
    non_d0_irred_bases := 0;
    d0_contact_bases := 0;
    contact_rank_hist := AssociativeArray();
    contact_type_hist := AssociativeArray();
    samples := [];

    stop := false;
    for rr in F do
        if stop then break; end if;
        if rr eq 0 then
            continue;
        end if;
        for tt in F do
            if stop then break; end if;
            if tt eq 0 or rr*tt eq 1 then
                continue;
            end if;
            for pp in F do
                base +:= 1;
                if MaxCurvesPerPrime gt 0 and base gt MaxCurvesPerPrime then
                    stop := true;
                    break;
                end if;
                f, g8, ellBase, dval := A8fFF(x, rr, pp, tt);
                if Degree(f) lt 5 or not IsSquarefree(f) or IsEvenSextic(f) then
                    continue;
                end if;
                smooth +:= 1;
                J := Jacobian(HyperellipticCurve(f));
                O := J!0;
                // Finite-field survey uses the unscaled A(8) model, so the
                // visible class is represented by -ellBase.  The older
                // rational searches multiply by the integer-model denominator
                // after scaling f by L^2.
                v8 := (-ellBase) mod g8;
                ok8 := true;
                try
                    D8 := elt<J | g8, v8>;
                    if 8*D8 ne O or 4*D8 eq O then
                        ok8 := false;
                    end if;
                catch e
                    ok8 := false;
                end try;
                if not ok8 then
                    continue;
                end if;
                order8 +:= 1;
                r2 := TwoRankFF(f);
                if r2 lt 2 then
                    continue;
                end if;
                wsplit +:= 1;
                if #J mod 3 ne 0 then
                    continue;
                end if;
                wsplit_j3 +:= 1;
                crows := ContactClassesFF(J, f);
                if #crows eq 0 then
                    continue;
                end if;
                contact_bases +:= 1;
                if dval eq 0 then
                    d0_contact_bases +:= 1;
                else
                    non_d0_contact_bases +:= 1;
                end if;
                has_irred := false;
                for row in crows do
                    contact_classes +:= 1;
                    u := row[1];
                    split := row[2];
                    if split then
                        split_classes +:= 1;
                        Inc(~contact_type_hist, "split");
                    else
                        irreducible_classes +:= 1;
                        Inc(~contact_type_hist, "irreducible");
                        has_irred := true;
                    end if;
                    AddSample(~samples, <rr, pp, tt, dval, r2, u, split>, SampleLimit);
                end for;
                Inc(~contact_rank_hist, Sprint(r2));
                if dval ne 0 and has_irred then
                    non_d0_irred_bases +:= 1;
                end if;
            end for;
        end for;
    end for;

    printf "%o %o %o %o %o %o %o %o %o %o %o %o %o\n",
        ell, base, smooth, order8, wsplit, wsplit_j3,
        contact_bases, contact_classes, split_classes, irreducible_classes,
        non_d0_contact_bases, non_d0_irred_bases, d0_contact_bases;
    print " rank_hist", Sort([ <k, contact_rank_hist[k]> : k in Keys(contact_rank_hist) ]);
    print " type_hist", Sort([ <k, contact_type_hist[k]> : k in Keys(contact_type_hist) ]);
    print " samples";
    for s in samples do
        print "  ", s;
    end for;
end for;

print "GENERAL_NON_D0_CONTACT_RECON_DONE";
quit;
