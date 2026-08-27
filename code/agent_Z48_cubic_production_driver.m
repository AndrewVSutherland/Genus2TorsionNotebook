//////////////////////////////////////////////////////////////////////
//  Agent Z48 cubic-contact production driver.
//
//  Production wrapper copied from agent_Z48_cubic_contact_route.m on
//  2026-07-02.  The only functional difference is robust initialization
//  of optional FixedR/FixedT parameters, so non-fixed partition scouts run
//  cleanly in this Magma invocation.
//
//  This is a bounded continuation after the RTHeight=4 A16 point-count
//  sweep.  It reuses the sign-aware A(8)->A(16) candidate equations, but
//  sends survivors through a genuine 3-part diagnostic:
//
//      h(x)^2 - f(x) = Lambda*(x^2 + U*x + V)^3.
//
//  For a smooth genus-2 curve with a verified rational order-16 point,
//  the Z/48 condition is equivalent to J(Q)[3] != 0.  The script therefore
//  uses exact TorsionSubgroup as the decisive 3-part test, and also runs a
//  bounded rational cubic-contact witness search as an explicit diagnostic.
//
//  Smoke:
//      magma -b Mode:=smoke ContactHeight:=2 \
//          code/agent_Z48_cubic_contact_route.m
//
//  Small bounded search:
//      magma -b Mode:=search RTHeight:=4 \
//          ExcludeRTHeight:=3 SearchBound:=2 PrimeBound:=19 MinGood:=1 \
//          MaxSlices:=8 ContactHeight:=1 \
//          code/agent_Z48_cubic_contact_route.m
//
//  Fixed known survivor slice:
//      magma -b Mode:=search FixedRNum:=-1 \
//          FixedRDen:=4 FixedTNum:=-1 FixedTDen:=4 SearchBound:=8 \
//          PrimeBound:=43 MinGood:=3 ContactHeight:=1 \
//          code/agent_Z48_cubic_contact_route.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Qq := Rationals();
Z := Integers();
Px<x> := PolynomialRing(Qq);
Pz<xz> := PolynomialRing(Z);

function StringParam(v)
    if Type(v) eq MonStgElt then
        return v;
    end if;
    return Sprint(v);
end function;

function BoolString(v)
    s := StringParam(v);
    return s in {"true", "True", "TRUE", "1", "yes", "Yes", "YES"};
end function;

if assigned Mode then
    Mode := StringParam(Mode);
else
    Mode := "smoke";
end if;

if assigned RTHeight and Type(RTHeight) eq MonStgElt then
    RTHeight := StringToInteger(RTHeight);
end if;
if not assigned RTHeight then
    RTHeight := 4;
end if;

if assigned ExcludeRTHeight and Type(ExcludeRTHeight) eq MonStgElt then
    ExcludeRTHeight := StringToInteger(ExcludeRTHeight);
end if;
if not assigned ExcludeRTHeight then
    ExcludeRTHeight := 3;
end if;

if assigned SearchBound and Type(SearchBound) eq MonStgElt then
    SearchBound := StringToInteger(SearchBound);
end if;
if not assigned SearchBound then
    SearchBound := 2;
end if;

if assigned PrimeBound and Type(PrimeBound) eq MonStgElt then
    PrimeBound := StringToInteger(PrimeBound);
end if;
if not assigned PrimeBound then
    PrimeBound := 43;
end if;

if assigned MinGood and Type(MinGood) eq MonStgElt then
    MinGood := StringToInteger(MinGood);
end if;
if not assigned MinGood then
    MinGood := 3;
end if;

if assigned GateMode then
    GateMode := StringParam(GateMode);
else
    GateMode := "48";
end if;

if assigned SliceMod and Type(SliceMod) eq MonStgElt then
    SliceMod := StringToInteger(SliceMod);
end if;
if not assigned SliceMod then
    SliceMod := 1;
end if;
if SliceMod lt 1 then
    error "SliceMod must be positive";
end if;

if assigned SliceClass and Type(SliceClass) eq MonStgElt then
    SliceClass := StringToInteger(SliceClass);
end if;
if not assigned SliceClass then
    SliceClass := 0;
end if;
if SliceClass lt 0 or SliceClass ge SliceMod then
    error "SliceClass must satisfy 0 <= SliceClass < SliceMod";
end if;

if assigned MaxSlices and Type(MaxSlices) eq MonStgElt then
    MaxSlices := StringToInteger(MaxSlices);
end if;
if not assigned MaxSlices then
    MaxSlices := 8;
end if;

if assigned MaxA16Roots and Type(MaxA16Roots) eq MonStgElt then
    MaxA16Roots := StringToInteger(MaxA16Roots);
end if;
if not assigned MaxA16Roots then
    MaxA16Roots := 20;
end if;

if assigned MaxExact and Type(MaxExact) eq MonStgElt then
    MaxExact := StringToInteger(MaxExact);
end if;
if not assigned MaxExact then
    MaxExact := 20;
end if;

if assigned ContactHeight and Type(ContactHeight) eq MonStgElt then
    ContactHeight := StringToInteger(ContactHeight);
end if;
if not assigned ContactHeight then
    ContactHeight := 2;
end if;

if assigned DoExact3 then
    DoExact3 := BoolString(DoExact3);
else
    DoExact3 := true;
end if;

if assigned ProgressSlices and Type(ProgressSlices) eq MonStgElt then
    ProgressSlices := StringToInteger(ProgressSlices);
end if;
if not assigned ProgressSlices then
    ProgressSlices := 1;
end if;

HaveFixedRNum := assigned FixedRNum;
HaveFixedRDen := assigned FixedRDen;
HaveFixedTNum := assigned FixedTNum;
HaveFixedTDen := assigned FixedTDen;

if not HaveFixedRNum then
    FixedRNum := 0;
end if;
if not HaveFixedRDen then
    FixedRDen := 1;
end if;
if not HaveFixedTNum then
    FixedTNum := 0;
end if;
if not HaveFixedTDen then
    FixedTDen := 1;
end if;

if HaveFixedRNum and Type(FixedRNum) eq MonStgElt then
    FixedRNum := StringToInteger(FixedRNum);
end if;
if HaveFixedRDen and Type(FixedRDen) eq MonStgElt then
    FixedRDen := StringToInteger(FixedRDen);
end if;
if HaveFixedTNum and Type(FixedTNum) eq MonStgElt then
    FixedTNum := StringToInteger(FixedTNum);
end if;
if HaveFixedTDen and Type(FixedTDen) eq MonStgElt then
    FixedTDen := StringToInteger(FixedTDen);
end if;
UseFixedSlice := HaveFixedRNum and HaveFixedRDen
                 and HaveFixedTNum and HaveFixedTDen;
if UseFixedSlice then
    if FixedRDen eq 0 or FixedTDen eq 0 then
        error "Fixed denominators must be nonzero";
    end if;
    FixedR := Qq!FixedRNum/Qq!FixedRDen;
    FixedT := Qq!FixedTNum/Qq!FixedTDen;
else
    FixedR := Qq!0;
    FixedT := Qq!0;
end if;

if assigned MaxGatePrint and Type(MaxGatePrint) eq MonStgElt then
    MaxGatePrint := StringToInteger(MaxGatePrint);
end if;
if not assigned MaxGatePrint then
    MaxGatePrint := 20;
end if;

procedure Bump(~A, key)
    if IsDefined(A, key) then
        A[key] +:= 1;
    else
        A[key] := 1;
    end if;
end procedure;

function HeightRationals(H)
    vals := [Qq!0];
    for den in [1..H] do
        for num in [-H..H] do
            if GCD(num, den) eq 1 then
                Append(~vals, Qq!num/Qq!den);
            end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

function CoeffSafe(f, i)
    if i gt Degree(f) then
        return Qq!0;
    end if;
    return Qq!Coefficient(f, i);
end function;

function NormalizeInvariants(inv)
    return Sort([Z!n : n in inv | Z!n ne 1]);
end function;

function TorsionOrder(inv)
    ninv := NormalizeInvariants(inv);
    if #ninv eq 0 then
        return 1;
    end if;
    return &*ninv;
end function;

function TorsionExponent(inv)
    ninv := NormalizeInvariants(inv);
    if #ninv eq 0 then
        return 1;
    end if;
    e := 1;
    for n in ninv do
        e := LCM(e, Z!n);
    end for;
    return e;
end function;

function HasThreePart(inv)
    return TorsionExponent(inv) mod 3 eq 0;
end function;

function HasOrder48(inv)
    return TorsionExponent(inv) mod 48 eq 0;
end function;

function SquareIntegralPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Pz!(L^2*f), L^2;
end function;

function A8Data(rv, pv, tv)
    r := Qq!rv;
    p := Qq!pv;
    t := Qq!tv;

    e := t^2 - 2*p*t/r;
    s := p - r^2;
    d := e + 2*p - r^2;
    lambda := r/t;
    u := p + r*t - 2*r;
    v := e + r^2 - r*p - r^2*t + 3*p*t - r*t^2;

    a := r^2 - lambda;
    b := 2*r*p - 2*lambda*(p + r*t) + 2*r*lambda;
    c := p^2 + 2*p*r^2 - r^4 - r^3*t - r*p^2/t
         - lambda*(r^2 + e)
         + 2*lambda*(r*p + r^2*t - 3*p*t + r*t^2);

    Qpoly := x^2 + d;
    q := a*x^2 + b*x + c;
    L := r*x + s;
    g8 := x^2 + u*x + v;
    f := q*(Qpoly^2 + q);
    ell8base := -(q + Qpoly*L);
    return f, g8, ell8base;
end function;

function A16Data(RVal, TVal, pv, muv, yv, Nq, zq)
    fQ, g8Q, ellBaseQ := A8Data(RVal, pv, TVal);
    ellQ := ellBaseQ + g8Q*(muv*x + Nq);
    WQ := x^2 + yv*x + zq;
    scale := muv^2 - 2*RVal*muv + RVal/TVal;
    return fQ, g8Q, ellBaseQ, ellQ, WQ, scale;
end function;

function VerifyA16Candidate(RVal, TVal, pv, muv, yv, Nq, zq)
    fQ, g8Q, ellBaseQ, ellQ, WQ, scale :=
        A16Data(RVal, TVal, pv, muv, yv, Nq, zq);
    if (ellQ^2 - fQ) mod g8Q ne 0 then
        return false, 0, false, fQ, "not_divisible";
    end if;
    SQ := ExactQuotient(ellQ^2 - fQ, g8Q);
    if SQ ne scale*WQ^2 then
        return false, 0, false, fQ, "not_square";
    end if;

    try
        Cc := HyperellipticCurve(fQ);
        J := Jacobian(Cc);
        ZJ := J!0;
        D8 := J![g8Q, -ellBaseQ mod g8Q];
        D16minus := J![WQ, -ellQ mod WQ];
        D16plus := J![WQ, ellQ mod WQ];

        if 8*D8 ne ZJ or 4*D8 eq ZJ then
            return false, 0, true, fQ, "d8_reject";
        end if;
        if 2*D16minus eq D8 and 16*D16minus eq ZJ and 8*D16minus ne ZJ then
            return true, -1, true, fQ, "certified";
        end if;
        if 2*D16plus eq D8 and 16*D16plus eq ZJ and 8*D16plus ne ZJ then
            return true, 1, true, fQ, "certified";
        end if;
    catch err
        return false, 0, true, fQ, "jacobian_error";
    end try;

    return false, 0, true, fQ, "d16_reject";
end function;

function GoodReductionPolynomial(f, ell)
    F := GF(ell);
    PF<xp> := PolynomialRing(F);
    try
        fp := PF![F!Coefficient(f, i) : i in [0..Degree(f)]];
    catch err
        return false, PF!0;
    end try;
    if Degree(fp) lt 5 or Discriminant(fp) eq 0 then
        return false, fp;
    end if;
    return true, fp;
end function;

function PointCountGate(f, primes, minGood, gateMode)
    used := [];
    gcdN := 0;
    good := 0;
    gateMod := (gateMode eq "3") select 3 else 48;

    for ell in primes do
        ok, fp := GoodReductionPolynomial(f, ell);
        if not ok then
            continue;
        end if;

        C := HyperellipticCurve(fp);
        Np := Z!#Jacobian(C);
        if gcdN eq 0 then
            gcdN := Np;
        else
            gcdN := GCD(gcdN, Np);
        end if;
        good +:= 1;
        Append(~used, <ell, Np, Np mod gateMod, gcdN>);

        if Np mod gateMod ne 0 then
            return false, "killed", ell, Np, good, gcdN, used;
        end if;
        if good ge minGood then
            return true, "passed", 0, 0, good, gcdN, used;
        end if;
    end for;

    return false, "insufficient_good", 0, 0, good, gcdN, used;
end function;

function EvalFraction(fr, val)
    num := Numerator(fr);
    den := Denominator(fr);
    denv := Evaluate(den, val);
    if denv eq 0 then
        return false, Qq!0;
    end if;
    return true, Qq!Evaluate(num, val)/Qq!denv;
end function;

function CubicContactForQ(f, U, V)
    q := x^2 + U*x + V;
    if Discriminant(q) eq 0 then
        return false, "disc_q_zero", "";
    end if;
    if Degree(GCD(q, f)) gt 0 then
        return false, "gcd_q_f", "";
    end if;

    q3 := q^3;
    a := [CoeffSafe(f, i) : i in [0..6]];
    c := [CoeffSafe(q3, i) : i in [0..6]];

    R<M> := PolynomialRing(Qq);
    FR := FieldOfFractions(R);
    H := FR!M;
    lam := H^2 - a[7];

    // Generic branch H != 0.  This is the branch relevant for the bounded
    // route; exact TorsionSubgroup remains the decisive 3-part test.
    H2 := (a[6] + lam*c[6])/(2*H);
    H1 := (a[5] + lam*c[5] - H2^2)/(2*H);
    H0 := (a[4] + lam*c[4] - 2*H2*H1)/(2*H);

    E2 := H1^2 + 2*H2*H0 - a[3] - lam*c[3];
    E1 := 2*H1*H0 - a[2] - lam*c[2];
    E0 := H0^2 - a[1] - lam*c[1];

    polys := [R!Numerator(E2), R!Numerator(E1), R!Numerator(E0)];
    G := R!0;
    for pol in polys do
        if pol ne 0 then
            if G eq 0 then
                G := pol;
            else
                G := GCD(G, pol);
            end if;
        end if;
    end for;

    if G eq 0 then
        return false, "positive_dim_generic_branch", "";
    end if;
    if Degree(G) le 0 then
        return false, "no_h3_roots", "";
    end if;

    for rt in Roots(G) do
        H3v := rt[1];
        if H3v eq 0 then
            continue;
        end if;
        ok2, H2v := EvalFraction(H2, H3v);
        ok1, H1v := EvalFraction(H1, H3v);
        ok0, H0v := EvalFraction(H0, H3v);
        if not (ok2 and ok1 and ok0) then
            continue;
        end if;
        lamv := H3v^2 - a[7];
        if lamv eq 0 then
            continue;
        end if;
        h := H3v*x^3 + H2v*x^2 + H1v*x + H0v;
        if h^2 - f eq lamv*q^3 then
            witness := Sprintf("U=%o V=%o H3=%o H2=%o H1=%o H0=%o Lambda=%o q=%o h=%o",
                U, V, H3v, H2v, H1v, H0v, lamv, q, h);
            return true, "contact", witness;
        end if;
    end for;

    return false, "roots_fail_verification", "";
end function;

function CubicContactHeightSearch(f, height)
    vals := HeightRationals(height);
    checked := 0;
    skipped := AssociativeArray();
    lastReason := "none";

    for U in vals do
        for V in vals do
            checked +:= 1;
            found, reason, witness := CubicContactForQ(f, U, V);
            if found then
                return true, checked, skipped, witness;
            end if;
            lastReason := reason;
            Bump(~skipped, reason);
        end for;
    end for;

    return false, checked, skipped, lastReason;
end function;

function ExactThreePart(f)
    try
        fInt, sqScale := SquareIntegralPolynomial(f);
        Ato, mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fInt)));
        inv := NormalizeInvariants(Invariants(Ato));
        return true, HasThreePart(inv), inv, TorsionOrder(inv), TorsionExponent(inv),
            sqScale, "";
    catch err
        return false, false, [], 0, 0, 0, Sprint(err`Object);
    end try;
end function;

function ThreePartDiagnostic(label, f, contactHeight, doExact)
    contactHit, contactChecked, contactSkipped, contactWitness :=
        CubicContactHeightSearch(f, contactHeight);
    exactOK := false;
    has3 := false;
    inv := [];
    torsOrder := 0;
    torsExponent := 0;
    sqScale := 0;
    exactReason := "not_run";

    if doExact then
        exactOK, has3, inv, torsOrder, torsExponent, sqScale, exactReason :=
            ExactThreePart(f);
    end if;

    printf "THREE_PART_DIAGNOSTIC label=%o contactHeight=%o contactChecked=%o contactHit=%o",
        label, contactHeight, contactChecked, contactHit;
    if contactHit then
        printf " contactWitness=\"%o\"", contactWitness;
    end if;
    if doExact then
        printf " exactOK=%o has3=%o torsion=%o order=%o exponent=%o squareScale=%o",
            exactOK, has3, inv, torsOrder, torsExponent, sqScale;
        if not exactOK then
            printf " exactReason=\"%o\"", exactReason;
        end if;
    end if;
    printf "\n";

    return (doExact select (exactOK and has3) else contactHit), contactHit,
        exactOK, has3, inv, torsOrder, torsExponent;
end function;

function SliceCategory(rOld, tOld)
    if rOld and tOld then
        return "completed";
    elif rOld then
        return "Tnew";
    elif tOld then
        return "Rnew";
    else
        return "bothnew";
    end if;
end function;

procedure PrintAssociativeCounts(A)
    for k in Sort([kk : kk in Keys(A)]) do
        printf "  %o : %o\n", k, A[k];
    end for;
end procedure;

TupleFormat := recformat<label, r, t, p, mu, y, N, z>;

smokeTuples := [
    rec<TupleFormat | label := "slice0_RT4_gate_survivor_false_positive",
        r := Qq!-1/4, t := Qq!-1/4, p := Qq!-41/144,
        mu := Qq!-1/2, y := Qq!-5/8, N := Qq!5/8, z := Qq!125/96>,
    rec<TupleFormat | label := "simple_rt_3_1over3_p2_mu9",
        r := Qq!3, t := Qq!1/3, p := Qq!2,
        mu := Qq!9, y := Qq!-4, N := Qq!-12, z := Qq!8/3>,
    rec<TupleFormat | label := "simple_rt_3_1over3_p2_mu27over11",
        r := Qq!3, t := Qq!1/3, p := Qq!2,
        mu := Qq!27/11, y := Qq!-8, N := Qq!-60/11, z := Qq!32/3>,
    rec<TupleFormat | label := "simple_rt_3_1over3_p34over9",
        r := Qq!3, t := Qq!1/3, p := Qq!34/9,
        mu := Qq!9/17, y := Qq!-4/21, N := Qq!-12/17,
        z := Qq!-1544/567>,
    rec<TupleFormat | label := "simple_rt_minus1_1over3_p1over3",
        r := Qq!-1, t := Qq!1/3, p := Qq!1/3,
        mu := Qq!3/7, y := Qq!4/3, N := Qq!8/7,
        z := Qq!8/3>,
    rec<TupleFormat | label := "simple_rt_minus1_1over2_pminus35over6",
        r := Qq!-1, t := Qq!1/2, p := Qq!-35/6,
        mu := Qq!1, y := Qq!-7, N := Qq!-13/2,
        z := Qq!77/4>,
    rec<TupleFormat | label := "simple_rt_3_1over2_p17over6",
        r := Qq!3, t := Qq!1/2, p := Qq!17/6,
        mu := Qq!1, y := Qq!-1/3, N := Qq!-7/6,
        z := Qq!-227/36>
];

procedure RunSmoke()
    primeList := [ell : ell in PrimesUpTo(PrimeBound) | ell notin {2,3}];

    print "Z48_CUBIC_CONTACT_SMOKE";
    printf "PrimeBound=%o MinGood=%o GateMode=%o ContactHeight=%o DoExact3=%o\n",
        PrimeBound, MinGood, GateMode, ContactHeight, DoExact3;

    rejected3 := 0;
    contactHits := 0;
    verified16 := 0;
    order48 := 0;

    for tup in smokeTuples do
        print "";
        printf "SMOKE_CANDIDATE label=%o r=%o t=%o p=%o mu=%o y=%o N=%o z=%o\n",
            tup`label, tup`r, tup`t, tup`p, tup`mu, tup`y, tup`N, tup`z;

        ok16, sign, squareOK, fQ, reason :=
            VerifyA16Candidate(tup`r, tup`t, tup`p, tup`mu, tup`y, tup`N, tup`z);
        if ok16 then
            verified16 +:= 1;
        end if;
        printf "A16_VERIFY label=%o squareOK=%o ok16=%o sign=%o reason=%o degree=%o disc_zero=%o factor_degrees=%o\n",
            tup`label, squareOK, ok16, sign, reason, Degree(fQ),
            Discriminant(fQ) eq 0,
            [<Degree(ff[1]), ff[2]> : ff in Factorization(fQ)];

        passGate, gateReason, killp, killN, goodCount, gcdN, used :=
            PointCountGate(fQ, primeList, MinGood, GateMode);
        printf "POINT_GATE label=%o pass=%o reason=%o killp=%o killN=%o good=%o gcd=%o used=%o\n",
            tup`label, passGate, gateReason, killp, killN, goodCount, gcdN, used;

        threePass, contactHit, exactOK, has3, inv, torsOrder, torsExponent :=
            ThreePartDiagnostic(tup`label, fQ, ContactHeight, DoExact3);
        if contactHit then
            contactHits +:= 1;
        end if;
        if not threePass then
            rejected3 +:= 1;
        end if;
        if ok16 and has3 then
            order48 +:= 1;
            printf "SMOKE_Z48_UNEXPECTED label=%o torsion=%o exponent=%o\n",
                tup`label, inv, torsExponent;
        end if;
    end for;

    printf "SMOKE_DONE candidates=%o verified16=%o rejectedByExact3=%o contactHeightHits=%o z48Unexpected=%o\n",
        #smokeTuples, verified16, rejected3, contactHits, order48;
end procedure;

procedure RunSearch()
    primeList := [ell : ell in PrimesUpTo(PrimeBound) | ell notin {2,3}];
    rvals := UseFixedSlice select [FixedR] else HeightRationals(RTHeight);
    tvals := UseFixedSlice select [FixedT] else HeightRationals(RTHeight);
    oldVals := (ExcludeRTHeight ge 0) select Seqset(HeightRationals(ExcludeRTHeight)) else {};
    searchVals := HeightRationals(SearchBound);

    print "Z48_CUBIC_CONTACT_SEARCH";
    printf "RTHeight=%o ExcludeRTHeight=%o SearchBound=%o PrimeBound=%o MinGood=%o GateMode=%o SliceMod=%o SliceClass=%o MaxSlices=%o MaxA16Roots=%o MaxExact=%o ContactHeight=%o DoExact3=%o UseFixedSlice=%o FixedR=%o FixedT=%o\n",
        RTHeight, ExcludeRTHeight, SearchBound, PrimeBound, MinGood,
        GateMode, SliceMod, SliceClass, MaxSlices, MaxA16Roots, MaxExact,
        ContactHeight, DoExact3, UseFixedSlice, FixedR, FixedT;
    print "primes", primeList;

    rawSlices := 0;
    completedSkipped := 0;
    eligibleSlices := 0;
    partitionSkipped := 0;
    maxSliceSkipped := 0;
    runSlices := 0;
    sliceBuildFail := 0;
    tested := 0;
    commonRootPairs := 0;
    rationalRoots := 0;
    singularReject := 0;
    nonsingularRoots := 0;
    pointGateReject := 0;
    pointGatePass := 0;
    threeChecked := 0;
    threeReject := 0;
    threePass := 0;
    contactHits := 0;
    exactTried := 0;
    squareReject := 0;
    d8Reject := 0;
    d16Reject := 0;
    certified := 0;
    z48Hits := 0;
    printedGate := 0;
    rootCapSkipped := 0;

    sliceCategoryHist := AssociativeArray();
    rootCategoryHist := AssociativeArray();
    nonsingularCategoryHist := AssociativeArray();
    firstKills := AssociativeArray();
    firstKillsByCategory := AssociativeArray();
    exact3TorsionHist := AssociativeArray();

    for RVal in rvals do
        if RVal eq 0 or RVal eq 1 then
            continue;
        end if;
        for TVal in tvals do
            if TVal eq 0 then
                continue;
            end if;

            rawSlices +:= 1;
            rOld := RVal in oldVals;
            tOld := TVal in oldVals;
            category := SliceCategory(rOld, tOld);

            if category eq "completed" then
                completedSkipped +:= 1;
                continue;
            end if;

            eligibleSlices +:= 1;
            if (eligibleSlices - 1) mod SliceMod ne SliceClass then
                partitionSkipped +:= 1;
                continue;
            end if;
            if MaxSlices gt 0 and runSlices ge MaxSlices then
                maxSliceSkipped +:= 1;
                continue;
            end if;

            runSlices +:= 1;
            Bump(~sliceCategoryHist, category);
            if ProgressSlices gt 0 and runSlices mod ProgressSlices eq 0 then
                printf "PROGRESS runSlices=%o rawSlices=%o tested=%o roots=%o nonsingular=%o gatePass=%o threePass=%o exact=%o certified=%o z48=%o\n",
                    runSlices, rawSlices, tested, rationalRoots,
                    nonsingularRoots, pointGatePass, threePass, exactTried,
                    certified, z48Hits;
            end if;

            try
                K<mu,y> := RationalFunctionField(Qq, 2);
                A<p,N,z> := PolynomialRing(K, 3);
                AX<X> := PolynomialRing(A);

                rv := K!RVal;
                tv := K!TVal;

                e := tv^2 - 2*p*tv/rv;
                s := p - rv^2;
                d := e + 2*p - rv^2;
                lambda := rv/tv;
                u := p + rv*tv - 2*rv;
                v := e + rv^2 - rv*p - rv^2*tv + 3*p*tv - rv*tv^2;

                a := rv^2 - lambda;
                b := 2*rv*p - 2*lambda*(p + rv*tv) + 2*rv*lambda;
                c := p^2 + 2*p*rv^2 - rv^4 - rv^3*tv - rv*p^2/tv
                     - lambda*(rv^2 + e)
                     + 2*lambda*(rv*p + rv^2*tv - 3*p*tv + rv*tv^2);

                Qpoly := X^2 + d;
                q := a*X^2 + b*X + c;
                L := rv*X + s;
                g8 := X^2 + u*X + v;
                f := q*(Qpoly^2 + q);
                ellBase := -(q + Qpoly*L);
                ell := ellBase + g8*(mu*X + N);
                S := ExactQuotient(ell^2 - f, g8);
                Clead := Coefficient(S, 4);
                W := X^2 + y*X + z;
                Diff := S - Clead*W^2;

                eq3 := Coefficient(Diff, 3);
                eq2 := Coefficient(Diff, 2);
                eq1 := Coefficient(Diff, 1);
                eq0 := Coefficient(Diff, 0);

                Nsol := -Coefficient(eq3, N, 0)/Coefficient(eq3, N, 1);

                B<p2,z2> := PolynomialRing(K, 2);
                NsolB := B!Coefficient(Nsol, 1, 0)
                         + (B!Coefficient(Nsol, 1, 1))*p2;
                hN := hom<A -> B | p2, NsolB, z2>;
                E2 := hN(eq2);
                E1 := hN(eq1);
                E0 := hN(eq0);
                zsol := -Coefficient(E2, z2, 0)/Coefficient(E2, z2, 1);

                C<p3> := PolynomialRing(K);
                zsolC := C!Coefficient(zsol, 1, 0)
                         + (C!Coefficient(zsol, 1, 1))*p3
                         + (C!Coefficient(zsol, 1, 2))*p3^2;
                NsolC := C!Coefficient(NsolB, 1, 0)
                         + (C!Coefficient(NsolB, 1, 1))*p3;
                hZ := hom<B -> C | p3, zsolC>;
                F1 := hZ(E1);
                F0 := hZ(E0);
            catch err
                sliceBuildFail +:= 1;
                printf "SLICE_BUILD_FAIL r=%o t=%o category=%o : %o\n",
                    RVal, TVal, category, err`Object;
                continue;
            end try;

            Qp<Pp> := PolynomialRing(Qq);
            for muv in searchVals do
                for yv in searchVals do
                    tested +:= 1;
                    if muv eq 0 then
                        continue;
                    end if;
                    if muv^2 - 2*RVal*muv + RVal/TVal eq 0 then
                        continue;
                    end if;

                    try
                        hK := hom<K -> Qq | muv, yv>;
                        hC := hom<C -> Qp | hK, Pp>;
                        F1q := hC(F1);
                        F0q := hC(F0);
                        Gp := GCD(F1q, F0q);
                    catch err
                        continue;
                    end try;

                    if Degree(Gp) le 0 then
                        continue;
                    end if;
                    commonRootPairs +:= 1;
                    roots := Roots(Gp);
                    for rt in roots do
                        if MaxA16Roots gt 0 and rationalRoots ge MaxA16Roots then
                            rootCapSkipped +:= 1;
                            continue;
                        end if;

                        pv := rt[1];
                        rationalRoots +:= 1;
                        Bump(~rootCategoryHist, category);
                        try
                            Nq := Evaluate(hC(NsolC), pv);
                            zq := Evaluate(hC(zsolC), pv);
                        catch err
                            continue;
                        end try;

                        fQ, g8Q, ellBaseQ := A8Data(RVal, pv, TVal);
                        if Degree(fQ) lt 5 or Discriminant(fQ) eq 0 then
                            singularReject +:= 1;
                            continue;
                        end if;

                        nonsingularRoots +:= 1;
                        Bump(~nonsingularCategoryHist, category);
                        passGate, gateReason, killp, killN, goodCount, gcdN, used :=
                            PointCountGate(fQ, primeList, MinGood, GateMode);
                        if not passGate then
                            pointGateReject +:= 1;
                            key := (gateReason eq "killed") select IntegerToString(killp)
                                   else gateReason;
                            Bump(~firstKills, key);
                            Bump(~firstKillsByCategory, category cat ":" cat key);
                            continue;
                        end if;

                        pointGatePass +:= 1;
                        if printedGate lt MaxGatePrint then
                            printf "GATE_PASS_BEFORE_CUBIC_CONTACT r=%o t=%o category=%o mu=%o y=%o p=%o N=%o z=%o good=%o gcd=%o used=%o\n",
                                RVal, TVal, category, muv, yv, pv, Nq, zq,
                                goodCount, gcdN, used;
                            printedGate +:= 1;
                        end if;

                        threeChecked +:= 1;
                        label := Sprintf("search_r=%o_t=%o_mu=%o_y=%o_p=%o",
                            RVal, TVal, muv, yv, pv);
                        pass3, contactHit, exactOK, has3, inv, torsOrder, torsExponent :=
                            ThreePartDiagnostic(label, fQ, ContactHeight, DoExact3);
                        if contactHit then
                            contactHits +:= 1;
                        end if;
                        if exactOK then
                            Bump(~exact3TorsionHist, Sprint(inv));
                        end if;
                        if not pass3 then
                            threeReject +:= 1;
                            continue;
                        end if;
                        threePass +:= 1;

                        if exactTried ge MaxExact then
                            continue;
                        end if;
                        exactTried +:= 1;

                        ok16, sign, squareOK, fCheck, reason :=
                            VerifyA16Candidate(RVal, TVal, pv, muv, yv, Nq, zq);
                        if not squareOK then
                            squareReject +:= 1;
                            continue;
                        end if;
                        if reason eq "d8_reject" then
                            d8Reject +:= 1;
                            continue;
                        end if;
                        if not ok16 then
                            d16Reject +:= 1;
                            continue;
                        end if;

                        certified +:= 1;
                        printf "CERTIFIED_A16_WITH_THREE_PART r=%o t=%o category=%o mu=%o y=%o p=%o N=%o z=%o sign=%o torsion=%o exponent=%o gate_gcd=%o gate_used=%o\n",
                            RVal, TVal, category, muv, yv, pv, Nq, zq,
                            sign, inv, torsExponent, gcdN, used;
                        printf "f=%o\n", fQ;
                        printf "factor_degrees=%o\n",
                            [<Degree(ff[1]), ff[2]> : ff in Factorization(fQ)];

                        if HasOrder48(inv) or has3 then
                            z48Hits +:= 1;
                            printf "Z48_CUBIC_CONTACT_HIT r=%o t=%o category=%o mu=%o y=%o p=%o torsion=%o exponent=%o\n",
                                RVal, TVal, category, muv, yv, pv, inv,
                                torsExponent;
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;

    printf "SEARCH_DONE rawSlices=%o completedSkipped=%o eligibleSlices=%o partitionSkipped=%o maxSliceSkipped=%o runSlices=%o sliceBuildFail=%o tested=%o commonRootPairs=%o rationalRoots=%o rootCapSkipped=%o singular=%o nonsingular=%o pointGateReject=%o pointGatePass=%o threeChecked=%o threeReject=%o threePass=%o contactHits=%o exactTried=%o squareReject=%o d8Reject=%o d16Reject=%o certified=%o z48Hits=%o\n",
        rawSlices, completedSkipped, eligibleSlices, partitionSkipped,
        maxSliceSkipped, runSlices, sliceBuildFail, tested, commonRootPairs,
        rationalRoots, rootCapSkipped, singularReject, nonsingularRoots,
        pointGateReject, pointGatePass, threeChecked, threeReject, threePass,
        contactHits, exactTried, squareReject, d8Reject, d16Reject,
        certified, z48Hits;
    print "SLICE_CATEGORIES_RUN";
    PrintAssociativeCounts(sliceCategoryHist);
    print "ROOT_CATEGORIES";
    PrintAssociativeCounts(rootCategoryHist);
    print "NONSINGULAR_CATEGORIES";
    PrintAssociativeCounts(nonsingularCategoryHist);
    print "FIRST_POINT_GATE_KILLS";
    PrintAssociativeCounts(firstKills);
    print "FIRST_POINT_GATE_KILLS_BY_CATEGORY";
    PrintAssociativeCounts(firstKillsByCategory);
    print "EXACT3_TORSION_HISTOGRAM";
    PrintAssociativeCounts(exact3TorsionHist);
end procedure;

print "AGENT_Z48_CUBIC_PRODUCTION_DRIVER";
printf "Mode=%o\n", Mode;

if Mode eq "smoke" then
    RunSmoke();
elif Mode eq "search" then
    RunSearch();
elif Mode eq "both" then
    RunSmoke();
    RunSearch();
else
    error "Mode must be smoke, search, or both";
end if;

print "DONE";
