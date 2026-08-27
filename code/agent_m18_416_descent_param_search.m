
//////////////////////////////////////////////////////////////////////
//  Parametrized descent search for P_R halving in M_1(8,4).
//
//  From agent_m18_416_descent_conditions.m (validated exactly):
//    P_R in 2J(Q)  <=>  C1: -c4*R = sq,  C2: A(-R) = sq,
//                       S_A, S_B: second-stage squares in the quadratic
//                       components of Q[T]/F.
//  With W = w^2 and K = -R(R^2-1):
//    C1 & (w rational)  <=>  w^2 - K s^2 = 1  (Pell conic, trivial point)
//    ==> FREE parametrization  w = (m^2 + K)/(m^2 - K).
//  On it, C2 collapses (mod squares) to
//    G(R,m) := (1-R^2)*(R*(2R+1) - W*(R+2))/2 = SQUARE,  W = w^2.
//
//  This script sweeps rational (R,m), tests G = sq, then the exact
//  second-stage descent, then certifies on the Jacobian and computes
//  torsion + simplicity.  Any full pass is a rational order-16 point
//  in M_1(8,4) (>= [2,16]); with the first tangent cover, [4,16].
//
//  Usage: magma -b hR:=40 hM:=40 NParts:=6 Part:=0 \
//             agent_m18_416_descent_param_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned hR then hR := 30;
elif Type(hR) eq MonStgElt then hR := StringToInteger(hR); end if;
if not assigned hM then hM := 30;
elif Type(hM) eq MonStgElt then hM := StringToInteger(hM); end if;
if not assigned NParts then NParts := 1;
elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0;
elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned progress then progress := 200000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

P<xq> := PolynomialRing(Q);

function IsSquareQ(qv)
    qv := Q!qv;
    if qv lt 0 then return false; end if;
    if qv eq 0 then return false; end if;   // exclude degenerate zeros here
    return IsSquare(Numerator(qv)) and IsSquare(Denominator(qv));
end function;

function FamilyData(Rv, wv)
    tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
    c4v := Rv + 2 + 4*tv;
    Av := xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4;
    Bv := c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv);
    return xq*Av*Bv, Av, Bv, c4v;
end function;

// exact second-stage descent test (norm conditions assumed prechecked)
function SecondStagePass(Rv, wv)
    fq, Av, Bv, c4v := FamilyData(Rv, wv);
    if Degree(fq) ne 5 or Discriminant(fq) eq 0 then return false; end if;
    XR := -c4v*Rv;
    Atq := P![Q!co : co in Coefficients(c4v^2*Evaluate(Av, xq/c4v))];
    Btq := P![Q!co : co in Coefficients(c4v*Evaluate(Bv, xq/c4v))];
    for gq in [Atq, Btq] do
        dsc := Discriminant(gq);
        if dsc eq 0 then return false; end if;
        if IsSquareQ(dsc) then
            for rt in Roots(gq) do
                val := XR - rt[1];
                if val eq 0 then return false; end if;
                if not IsSquareQ(val) then return false; end if;
            end for;
        else
            K<th> := NumberField(gq);
            if not IsSquare(XR - th) then return false; end if;
        end if;
    end for;
    return true;
end function;

function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do v := Evaluate(fp, xx);
        if v eq 0 then cnt +:= 1; elif IsSquare(v) then cnt +:= 2; end if; end for;
    if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if;
    return cnt;
end function;

function SimplicityCertificate(fInt)
    RT := PolynomialRing(Q); T := RT.1;
    dsc := Discriminant(fInt);
    for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        if (Z!LeadingCoefficient(fInt)) mod pp eq 0 then continue; end if;
        if (Z!Numerator(dsc)) mod pp eq 0 then continue; end if;
        PF := PolynomialRing(GF(pp));
        fp := PF![GF(pp)!co : co in Coefficients(fInt)];
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        PF2 := PolynomialRing(GF(pp^2));
        fp2 := PF2![GF(pp^2)!co : co in Coefficients(fInt)];
        a1 := pp + 1 - CountCurve(fp);
        a2 := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;
        chi := T^4 - a1*T^3 + a2*T^2 - a1*pp*T + pp^2;
        if not IsIrreducible(chi) then continue; end if;
        K := NumberField(chi); pi := K.1; drop := false;
        for n in [2..12] do
            if Degree(MinimalPolynomial(pi^n)) lt 4 then drop := true; break; end if;
        end for;
        if not drop then return true, pp, chi; end if;
    end for;
    return false, 0, RT!0;
end function;

function Has416(invs)
    vals := Reverse(Sort([Valuation(n, 2) : n in invs]));
    return #vals ge 2 and vals[1] ge 4 and vals[2] ge 2;
end function;

function RationalParametersOfHeight(Bd)
    vals := [];
    for den in [1..Bd] do for num in [-Bd..Bd] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

printf "PARAM DESCENT SEARCH hR=%o hM=%o Part=%o/%o\n", hR, hM, Part, NParts;
Rparams := RationalParametersOfHeight(hR);
Mparams := RationalParametersOfHeight(hM);

tested := 0; wOK := 0; c2pass := 0; sstage := 0; halvable := 0; hits := 0;
torsHist := AssociativeArray();
seenW := {};
ridx := 0;

for Rv in Rparams do
    ridx +:= 1;
    if (ridx mod NParts) ne Part then continue; end if;
    if Rv in {Q!0, Q!1, Q!-1} then continue; end if;
    K := -2*Rv*(Rv^2 - 1);
    for mv in Mparams do
        tested +:= 1;
        if tested mod progress eq 0 then
            printf "PROGRESS tested=%o wOK=%o C2pass=%o sstage=%o halvable=%o hits=%o\n",
                tested, wOK, c2pass, sstage, halvable, hits;
        end if;
        den := mv^2 - K;
        if den eq 0 then continue; end if;
        wv := (mv^2 + K)/den;
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        if <Rv,wv> in seenW then continue; end if;
        Include(~seenW, <Rv,wv>);
        wOK +:= 1;
        Wv := wv^2;
        // C2 = A(-R) = -2*R*(R-1)^2*Qfac/(w^2-1); on the Pell
        // parametrization (W-1) == K = -2*R*(R^2-1) mod squares, so
        // C2 == 2*(R^2-1)*(R*(2R+1) - W*(R+2)) mod squares.
        G := 2*(Rv^2-1)*(Rv*(2*Rv+1) - Wv*(Rv+2));
        if not IsSquareQ(G) then continue; end if;
        c2pass +:= 1;
        // exact second-stage
        if not SecondStagePass(Rv, wv) then continue; end if;
        sstage +:= 1;
        // exact Jacobian certification
        fq, Av, Bv, c4v := FamilyData(Rv, wv);
        L := 1;
        for i in [0..Degree(fq)] do L := LCM(L, Denominator(Coefficient(fq, i))); end for;
        fI := P!(L^2*fq);
        J := Jacobian(HyperellipticCurve(fI));
        Qf := Rv^2 - (Q!1/2)*Rv*wv^2 + (Q!1/2)*Rv - wv^2;
        YRv := -2*Rv*(Rv-1)^2*Qf/(wv^2-1);
        if Evaluate(fI, -Rv) ne (L*YRv)^2 then continue; end if;
        PR := J![xq + Rv, P!(L*YRv)];
        ok, half := IsDivisibleBy(PR, 2);
        printf "DESCENT_FULL_PASS R=%o w=%o exact=%o\n", Rv, wv, ok;
        if not ok then
            printf "  WARNING criterion passed but IsDivisibleBy failed\n";
            continue;
        end if;
        halvable +:= 1;
        ordh := Order(half);
        invs := Invariants(TorsionSubgroup(J));
        key := Sprint(invs);
        if IsDefined(torsHist,key) then torsHist[key] +:= 1; else torsHist[key] := 1; end if;
        printf "HALVABLE R=%o w=%o half_order=%o torsion=%o\n", Rv, wv, ordh, invs;
        printf "  f = %o\n", fI;
        issimple, pp, chi := SimplicityCertificate(fI);
        printf "  simplicity: %o\n",
            issimple select Sprintf("SIMPLE at p=%o chi=%o", pp, chi)
                     else "no certificate in small primes";
        if Has416(invs) then
            hits +:= 1;
            printf "HIT_416 R=%o w=%o torsion=%o simple=%o\n", Rv, wv, invs, issimple;
        end if;
    end for;
end for;

printf "SEARCH_DONE tested=%o wOK=%o C2pass=%o sstage=%o halvable=%o hits416=%o\n",
    tested, wOK, c2pass, sstage, halvable, hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsHist)]) do
    printf "  %o : %o\n", key, torsHist[key];
end for;
print "DONE";
quit;
