//////////////////////////////////////////////////////////////////////
//  Agent Z/48 simultaneous A(16)+3 search.
//
//  This searches the sign-aware square-root presentation of the
//  A(8)->A(16) cover.  For each fixed (r,t), the x^3 and x^2 square
//  equations solve N and z; the remaining x^1 and x^0 equations are
//  polynomials F1(p), F0(p) over Q(mu,y).  Small rational (mu,y) are
//  enumerated and common rational p-roots give A(16) candidates.
//
//  The Z/48 simultaneity is imposed before expensive exact Jacobian
//  certification: a certified order-16 point plus rational 3-torsion
//  would force 48 | #J(F_l) at every good prime l != 3.  Equivalently
//  the 3-part gate alone is necessary, but the default gate is 48 since
//  candidates are meant to lie on A(16).
//
//  Typical runs:
//      magma -b RTHeight:=2 SearchBound:=6 MinGood:=2 \
//          code/agent_Z48_simultaneous_A16_plus3.m
//
//      magma -b RTHeight:=3 SearchBound:=10 PrimeBound:=43 MinGood:=3 \
//          code/agent_Z48_simultaneous_A16_plus3.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Qq := Rationals();
Z := Integers();
Px<x> := PolynomialRing(Qq);
Pz<xz> := PolynomialRing(Z);

if assigned RTHeight and Type(RTHeight) eq MonStgElt then
    RTHeight := StringToInteger(RTHeight);
end if;
if not assigned RTHeight then
    RTHeight := 2;
end if;

if assigned SearchBound and Type(SearchBound) eq MonStgElt then
    SearchBound := StringToInteger(SearchBound);
end if;
if not assigned SearchBound then
    SearchBound := 6;
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

if assigned GateMode and Type(GateMode) ne MonStgElt then
    GateMode := Sprint(GateMode);
end if;
if not assigned GateMode then
    GateMode := "48";
end if;

if assigned MaxGatePrint and Type(MaxGatePrint) eq MonStgElt then
    MaxGatePrint := StringToInteger(MaxGatePrint);
end if;
if not assigned MaxGatePrint then
    MaxGatePrint := 20;
end if;

if assigned MaxExact and Type(MaxExact) eq MonStgElt then
    MaxExact := StringToInteger(MaxExact);
end if;
if not assigned MaxExact then
    MaxExact := 200;
end if;

if assigned ProgressSlices and Type(ProgressSlices) eq MonStgElt then
    ProgressSlices := StringToInteger(ProgressSlices);
end if;
if not assigned ProgressSlices then
    ProgressSlices := 10;
end if;

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

function ContainsZ48(inv)
    ninv := NormalizeInvariants(inv);
    return &or[n mod 48 eq 0 : n in ninv];
end function;

function SquareIntegralPolynomial(f)
    L := 1;
    for c in Coefficients(f) do
        L := LCM(L, Denominator(c));
    end for;
    return Pz!(L^2*f), L^2;
end function;

function VerifyA16Candidate(RVal, TVal, pv, muv, yv, Nq, zq)
    fQ, g8Q, ellBaseQ := A8Data(RVal, pv, TVal);
    ellQ := ellBaseQ + g8Q*(muv*x + Nq);
    if (ellQ^2 - fQ) mod g8Q ne 0 then
        return false, 0, false, fQ, "not_divisible";
    end if;
    SQ := ExactQuotient(ellQ^2 - fQ, g8Q);
    WQ := x^2 + yv*x + zq;
    scale := muv^2 - 2*RVal*muv + RVal/TVal;
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

primeList := [ell : ell in PrimesUpTo(PrimeBound) | ell notin {2,3}];
rvals := HeightRationals(RTHeight);
tvals := HeightRationals(RTHeight);
searchVals := HeightRationals(SearchBound);

print "AGENT_Z48_SIMULTANEOUS_A16_PLUS3";
printf "RTHeight=%o SearchBound=%o PrimeBound=%o MinGood=%o GateMode=%o MaxExact=%o\n",
    RTHeight, SearchBound, PrimeBound, MinGood, GateMode, MaxExact;
print "primes", primeList;

sliceIndex := 0;
sliceBuildFail := 0;
tested := 0;
commonRootPairs := 0;
rationalRoots := 0;
singularReject := 0;
pointGateReject := 0;
pointGatePass := 0;
exactTried := 0;
squareReject := 0;
d8Reject := 0;
d16Reject := 0;
certified := 0;
z48Hits := 0;
printedGate := 0;
firstKills := AssociativeArray();
certifiedTorsionHist := AssociativeArray();

for RVal in rvals do
    if RVal eq 0 or RVal eq 1 then
        continue;
    end if;
    for TVal in tvals do
        if TVal eq 0 then
            continue;
        end if;

        sliceIndex +:= 1;
        if ProgressSlices gt 0 and sliceIndex mod ProgressSlices eq 0 then
            printf "PROGRESS slices=%o tested=%o roots=%o gatePass=%o exact=%o certified=%o z48=%o\n",
                sliceIndex, tested, rationalRoots, pointGatePass,
                exactTried, certified, z48Hits;
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
            printf "SLICE_BUILD_FAIL r=%o t=%o : %o\n", RVal, TVal, err`Object;
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
                    pv := rt[1];
                    rationalRoots +:= 1;
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

                    passGate, gateReason, killp, killN, goodCount, gcdN, used :=
                        PointCountGate(fQ, primeList, MinGood, GateMode);
                    if not passGate then
                        pointGateReject +:= 1;
                        key := (gateReason eq "killed") select IntegerToString(killp)
                               else gateReason;
                        if IsDefined(firstKills, key) then
                            firstKills[key] +:= 1;
                        else
                            firstKills[key] := 1;
                        end if;
                        continue;
                    end if;

                    pointGatePass +:= 1;
                    if printedGate lt MaxGatePrint then
                        printf "GATE_PASS_A16_PLUS3 r=%o t=%o mu=%o y=%o p=%o N=%o z=%o good=%o gcd=%o used=%o\n",
                            RVal, TVal, muv, yv, pv, Nq, zq,
                            goodCount, gcdN, used;
                        printedGate +:= 1;
                    end if;

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
                    printf "CERTIFIED_A16_AFTER_Z48_GATE r=%o t=%o mu=%o y=%o p=%o N=%o z=%o sign=%o gate_gcd=%o gate_used=%o\n",
                        RVal, TVal, muv, yv, pv, Nq, zq, sign, gcdN, used;
                    printf "f=%o\n", fQ;
                    printf "factor_degrees=%o\n",
                        [<Degree(ff[1]), ff[2]> : ff in Factorization(fQ)];

                    try
                        fInt, sqScale := SquareIntegralPolynomial(fQ);
                        Ato, mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fInt)));
                        inv := NormalizeInvariants(Invariants(Ato));
                        key := Sprint(inv);
                        if IsDefined(certifiedTorsionHist, key) then
                            certifiedTorsionHist[key] +:= 1;
                        else
                            certifiedTorsionHist[key] := 1;
                        end if;
                        printf "integral_square_scale=%o torsion_invariants=%o torsion_order=%o\n",
                            sqScale, inv, TorsionOrder(inv);
                        if ContainsZ48(inv) then
                            z48Hits +:= 1;
                            printf "Z48_HIT r=%o t=%o mu=%o y=%o p=%o torsion=%o fInt=%o\n",
                                RVal, TVal, muv, yv, pv, inv, fInt;
                        end if;
                    catch err
                        printf "torsion computation failed: %o\n", err`Object;
                    end try;
                end for;
            end for;
        end for;
    end for;
end for;

printf "SEARCH_DONE slices=%o sliceBuildFail=%o tested=%o commonRootPairs=%o rationalRoots=%o singular=%o pointGateReject=%o pointGatePass=%o exactTried=%o squareReject=%o d8Reject=%o d16Reject=%o certified=%o z48Hits=%o\n",
    sliceIndex, sliceBuildFail, tested, commonRootPairs, rationalRoots,
    singularReject, pointGateReject, pointGatePass, exactTried,
    squareReject, d8Reject, d16Reject, certified, z48Hits;
print "FIRST_POINT_GATE_KILLS";
for k in Sort([kk : kk in Keys(firstKills)]) do
    printf "  %o : %o\n", k, firstKills[k];
end for;
print "CERTIFIED_TORSION_HISTOGRAM";
for k in Sort([kk : kk in Keys(certifiedTorsionHist)]) do
    printf "  %o : %o\n", k, certifiedTorsionHist[k];
end for;
print "DONE";

