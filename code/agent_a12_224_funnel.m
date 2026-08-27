
//////////////////////////////////////////////////////////////////////
//  Z/2 x Z/24 second-stage descent funnel on the A(12) chart.
//
//  Locus: (p, z, v) with r = (v^2 - p^2)/p, which makes the necessary
//  norm conditions N1 == N2 == p(p+r) = v^2 hold identically
//  (agent_a12_224_descent_setup.m, validated).
//
//  Per point, D = P12 (and the translate P12 + T_R, T_R = [R_m, 0]):
//  delta(D) trivial in A*/(A*^2 Q*)  <=>  exists q in Q*:
//     w_R/q in K_R*^2  and  w_F/q in K_F*^2,
//  where w_R = (u4*u6 mod R_m)(theta_R), w_F = (u4*u6 mod F_m)(theta_F).
//  The K_R side determines q up to <squares, d_R>: candidates
//     q in {2(a+n), 2(a-n)} * {1, d_R},   n = sqrt(N_{K_R}(w_R)).
//  Only candidates passing the (cheap) K_R test reach the quartic
//  IsSquare test in K_F.  Funnel passes go to exact IsDivisibleBy +
//  TorsionSubgroup + simplicity certificate.
//
//  Modes:
//    Validate:=true  -- small box, compare funnel vs IsDivisibleBy
//                       (funnel must never reject a halvable point);
//    default         -- search box pH/zH/vH with NParts/Part.
//
//  Usage:
//    magma -b Validate:=true agent_a12_224_funnel.m
//    magma -b pH:=20 zH:=20 vH:=20 NParts:=6 Part:=0 agent_a12_224_funnel.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
// hard memory cap per process (machine safety); override with MemGB
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned pH then pH := 10; elif Type(pH) eq MonStgElt then pH := StringToInteger(pH); end if;
if not assigned zH then zH := 10; elif Type(zH) eq MonStgElt then zH := StringToInteger(zH); end if;
if not assigned vH then vH := 10; elif Type(vH) eq MonStgElt then vH := StringToInteger(vH); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned Validate then Validate := false;
elif Type(Validate) eq MonStgElt then Validate := Validate in {"true","True","1","yes"}; end if;
if not assigned progress then progress := 100000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false, 0; end if;
    okn, sn := IsSquare(Numerator(qv));
    okd, sd := IsSquare(Denominator(qv));
    if okn and okd then return true, sn/sd; end if;
    return false, 0;
end function;

function SqfreePartQ(qv)
    qv := Q!qv;
    if qv eq 0 then return Q!0; end if;
    n := Z!(Numerator(qv)*Denominator(qv));
    sgn := Sign(n); n := Abs(n);
    sf := 1;
    for pe in Factorization(n) do
        if IsOdd(pe[2]) then sf *:= pe[1]; end if;
    end for;
    return Q!(sgn*sf);
end function;

function ChartData(pv, zv, rv)
    if pv eq 0 or zv eq 0 then return false, _, _, _, _, _; end if;
    sv := (zv^2 - 4*pv^2 + 1)/(2*zv);
    if sv^2 eq 1 then return false, _, _, _, _, _; end if;
    tv := (zv^2 + 4*pv^2 - 1)^2/(8*pv^2*zv);
    muv := ((sv^2 - 1)*(2*pv*rv + 1) - pv^2*(2*sv*tv - 4))/(4*pv^3);
    lav := (4 - muv^2)*pv^2/(sv^2 - 1);
    if lav eq 0 then return false, _, _, _, _, _; end if;
    T1v := pv*x + rv;
    Rv := (T1v^2 + x - 1)/lav;
    ellv := sv*x + tv;
    Qv := 2*T1v + muv*Rv;
    Fv := Rv*x^2 + 4*(Rv + x - 1)*(Rv - 1);
    fv := Rv*Fv;
    if Degree(fv) ne 6 or Discriminant(fv) eq 0 then return false, _, _, _, _, _; end if;
    return true, fv, Rv, Qv, ellv, Fv;
end function;

// quadratic-component data for w = (u mod g), g = x^2 + g1 x + g0 monic:
// representing w = A + B*theta.  Returns a-part (A - B*g1/2), b-part B,
// norm N(w), disc d = g1^2 - 4 g0.
function QuadData(w, g)
    A := Coefficient(w, 0); B := Coefficient(w, 1);
    g1 := Coefficient(g, 1); g0 := Coefficient(g, 0);
    Nw := A^2 - A*B*g1 + B^2*g0;
    aw := A - B*g1/2;
    d := g1^2 - 4*g0;
    return aw, B, Nw, d;
end function;

// is w/q a square in the quadratic algebra Q[x]/g (g monic quadratic)?
// (assumes N(w/q) = N(w)/q^2 is a rational square)
function QuadSquareTest(w, g, q)
    aw, bw, Nw, d := QuadData(w, g);
    okn, n := IsSquareQ(Nw/q^2);
    if not okn then return false; end if;
    // (w/q) square <=> 2(aw/q + n) or 2(aw/q - n) in Q*^2  (d-twist folded
    // into the q candidates)
    for sgn in [1,-1] do
        val := 2*(aw/q + sgn*n);
        if val ne 0 then
            ok, _ := IsSquareQ(val);
            if ok then return true; end if;
        else
            // w/q = +- n means b-part 0: rational element: square iff n>0 sq
            ok, _ := IsSquareQ(aw/q);
            if ok then return true; end if;
        end if;
    end for;
    return false;
end function;

// full second-stage descent test for a class with monic Mumford-u product
// UD (deg <= 4 poly; the product u4*u6 reduced appropriately):
// returns true if exists q: (UD mod Rm)/q square in K_R and
// (UD mod Fm)/q square in K_F.
function DescentTrivial(UD, Rm, Fm)
    wR := UD mod Rm;
    if Degree(wR) lt 0 then return false; end if;   // gcd issue
    aw, bw, Nw, d := QuadData(wR, Rm);
    if Nw eq 0 then return false; end if;
    okn, n := IsSquareQ(Nw);
    if not okn then return false; end if;
    // q candidates (squarefree representatives)
    cands := {};
    for sgn in [1,-1] do
        val := 2*(aw + sgn*n);
        if val ne 0 then
            Include(~cands, SqfreePartQ(val));
            Include(~cands, SqfreePartQ(val*d));
        end if;
    end for;
    if bw eq 0 then
        // w rational: q = w itself (up to squares) and d-twist
        Include(~cands, SqfreePartQ(aw));
        Include(~cands, SqfreePartQ(aw*d));
    end if;
    if #cands eq 0 then return false; end if;
    // K_F data
    wF := UD mod Fm;
    if Degree(wF) lt 0 then return false; end if;
    KF := NumberField(Fm);
    thF := KF.1;
    wFel := Evaluate(wF, thF);
    if wFel eq 0 then return false; end if;
    for q in cands do
        if q eq 0 then continue; end if;
        if not QuadSquareTest(wR, Rm, q) then continue; end if;
        if IsSquare(wFel/q) then return true; end if;
    end for;
    return false;
end function;

// exact side helpers
function BuildJacobian(fv)
    L := 1;
    for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
    fI := P!(L^2*fv);
    return Jacobian(HyperellipticCurve(fI)), L;
end function;

function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
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
        for nn in [2..12] do
            if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
        end for;
        if not drop then return true, pp, chi; end if;
    end for;
    return false, 0, RT!0;
end function;

// process one (pv, zv, rv): returns status string
function ProcessPoint(pv, zv, rv, doExact)
    ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
    if not ok then return "degenerate", _; end if;
    Rm := Rv/LeadingCoefficient(Rv);
    Fm := Fv/LeadingCoefficient(Fv);
    u4 := Qv/LeadingCoefficient(Qv);
    u6 := (Rv + x - 1)/LeadingCoefficient(Rv + x - 1);
    UD := u4*u6;
    // gcd checks
    if GCD(UD, Rm*Fm) ne 1 then return "gcd", _; end if;
    pass := DescentTrivial(UD, Rm, Fm);
    // translate P12 + T_R: compute its Mumford u via Jacobian arithmetic
    passT := false;
    if not pass then
        J, L := BuildJacobian(fv);
        try
            P4 := J![u4, (L*Rv*ellv) mod u4];
            P6 := J![u6, (L*x*Rv) mod u6];
            TR := J![Rm, P!0];
            D := P4 + P6 + TR;
            uD := Eltseq(D)[1];
            uD := uD/LeadingCoefficient(uD);
            if GCD(uD, Rm*Fm) eq 1 then
                // necessary norms for the translate are not free: quick gates
                wR := uD mod Rm;
                _, _, NwR, _ := QuadData(wR, Rm);
                oknr, _ := IsSquareQ(NwR);
                if oknr then
                    passT := DescentTrivial(uD, Rm, Fm);
                end if;
            end if;
        catch e
            passT := false;
        end try;
    end if;
    if not (pass or passT) then return "reject", _; end if;
    if not doExact then return "funnel-pass", _; end if;

    // exact certification
    J, L := BuildJacobian(fv);
    O := J!0;
    P4 := J![u4, (L*Rv*ellv) mod u4];
    P6 := J![u6, (L*x*Rv) mod u6];
    P12 := P4 + P6;
    if not (12*P12 eq O and &and[nn*P12 ne O : nn in [1..11]]) then
        return "not-order-12", _;
    end if;
    TR := J![Rm, P!0];
    for D in [P12, P12 + TR] do
        if not (12*D eq O) then continue; end if;
        div2, half := IsDivisibleBy(D, 2);
        if div2 then
            invs := Invariants(TorsionSubgroup(J));
            return "HALVABLE", <invs, Order(half)>;
        end if;
    end for;
    return "funnel-only", _;
end function;

function HeightRationals(H)
    vals := [];
    for den in [1..H] do for num in [-H..H] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

if Validate then
    print "== VALIDATION: funnel vs exact IsDivisibleBy on a small box ==";
    nchk := 0; nhalv := 0; nreject_halv := 0; npass := 0;
    for pv in HeightRationals(4) do
        if pv eq 0 then continue; end if;
        for zv in HeightRationals(3) do
            if zv eq 0 then continue; end if;
            for rv in HeightRationals(3) do
                ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
                if not ok then continue; end if;
                Rm := Rv/LeadingCoefficient(Rv);
                Fm := Fv/LeadingCoefficient(Fv);
                u4 := Qv/LeadingCoefficient(Qv);
                u6 := (Rv + x - 1)/LeadingCoefficient(Rv + x - 1);
                if GCD(u4*u6, Rm*Fm) ne 1 then continue; end if;
                // exact
                J, L := BuildJacobian(fv);
                O := J!0;
                try
                    P4 := J![u4, (L*Rv*ellv) mod u4];
                    P6 := J![u6, (L*x*Rv) mod u6];
                catch e continue; end try;
                P12 := P4 + P6;
                if not (12*P12 eq O and &and[nn*P12 ne O : nn in [1..11]]) then continue; end if;
                div2, _ := IsDivisibleBy(P12, 2);
                nchk +:= 1;
                funnel := DescentTrivial(u4*u6, Rm, Fm);
                if funnel then npass +:= 1; end if;
                if div2 then
                    nhalv +:= 1;
                    if not funnel then
                        nreject_halv +:= 1;
                        printf "FUNNEL FALSE-REJECT p=%o z=%o r=%o\n", pv, zv, rv;
                    end if;
                end if;
            end for;
        end for;
    end for;
    printf "validation: checked=%o funnel-pass=%o exactly-halvable=%o FALSE-REJECTS=%o\n",
        nchk, npass, nhalv, nreject_halv;

    // positive controls: D = 2*P12 is halvable BY CONSTRUCTION on every
    // chart curve; the funnel must pass its Mumford class.
    print "== positive controls: funnel on 2*P12 (halvable by construction) ==";
    npos := 0; nposreject := 0;
    for pv in HeightRationals(3) do
        if pv eq 0 then continue; end if;
        for zv in HeightRationals(2) do
            if zv eq 0 then continue; end if;
            for rv in HeightRationals(2) do
                ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
                if not ok then continue; end if;
                Rm := Rv/LeadingCoefficient(Rv);
                Fm := Fv/LeadingCoefficient(Fv);
                u4 := Qv/LeadingCoefficient(Qv);
                u6 := (Rv + x - 1)/LeadingCoefficient(Rv + x - 1);
                J, L := BuildJacobian(fv);
                O := J!0;
                try
                    P4 := J![u4, (L*Rv*ellv) mod u4];
                    P6 := J![u6, (L*x*Rv) mod u6];
                catch e continue; end try;
                P12 := P4 + P6;
                if not (12*P12 eq O and &and[nn*P12 ne O : nn in [1..11]]) then continue; end if;
                D2 := 2*P12;    // order 6, halvable by P12
                uD := Eltseq(D2)[1];
                if Degree(uD) lt 1 then continue; end if;
                uD := uD/LeadingCoefficient(uD);
                if GCD(uD, Rm*Fm) ne 1 then continue; end if;
                npos +:= 1;
                if not DescentTrivial(uD, Rm, Fm) then
                    nposreject +:= 1;
                    printf "POSITIVE-CONTROL REJECT p=%o z=%o r=%o\n", pv, zv, rv;
                end if;
            end for;
        end for;
    end for;
    printf "positive controls: %o tested, %o false-rejects\n", npos, nposreject;
    print "DONE";
    quit;
end if;

// ---- search over the parametrized conic locus ----
printf "SEARCH pH=%o zH=%o vH=%o Part=%o/%o  (r = (v^2-p^2)/p)\n",
    pH, zH, vH, Part, NParts;
pvals := HeightRationals(pH);
zvals := HeightRationals(zH);
vvals := HeightRationals(vH);
tested := 0; funnelpass := 0; halvable := 0;
pidx := 0;
for pv in pvals do
    pidx +:= 1;
    if (pidx mod NParts) ne Part then continue; end if;
    if pv eq 0 then continue; end if;
    for vv in vvals do
        rv := (vv^2 - pv^2)/pv;
        for zv in zvals do
            if zv eq 0 then continue; end if;
            tested +:= 1;
            if tested mod progress eq 0 then
                printf "PROGRESS tested=%o funnelpass=%o halvable=%o\n",
                    tested, funnelpass, halvable;
            end if;
            status, data := ProcessPoint(pv, zv, rv, true);
            if status eq "HALVABLE" then
                halvable +:= 1;
                printf "HALVABLE p=%o z=%o r=%o v=%o torsion=%o half_order=%o\n",
                    pv, zv, rv, vv, data[1], data[2];
                ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
                L := 1;
                for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
                fI := P!(L^2*fv);
                printf "  f = %o\n", fI;
                issimple, pp, chi := SimplicityCertificate(fI);
                printf "  simplicity: %o\n",
                    issimple select Sprintf("SIMPLE p=%o chi=%o", pp, chi)
                             else "none in small primes";
                inv := data[1];
                vals2 := Reverse(Sort([Valuation(nn,2) : nn in inv]));
                if #inv ge 2 and vals2[1] ge 3 and inv[#inv] mod 24 eq 0 then
                    printf "TARGET_2_24 p=%o z=%o r=%o torsion=%o simple=%o\n",
                        pv, zv, rv, inv, issimple;
                end if;
            elif status eq "funnel-pass" or status eq "funnel-only" then
                funnelpass +:= 1;
            end if;
        end for;
    end for;
end for;
printf "SEARCH_DONE tested=%o funnelpass=%o halvable=%o\n", tested, funnelpass, halvable;
print "DONE";
quit;
