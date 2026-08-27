//////////////////////////////////////////////////////////////////////
// Mordell-Weil sieve for the remaining R = -8 ELS fiber.
//
// The corrected S_B V_4 model has quotient
//
//     E_m : a^2 = -14*l^3 - (38/3)*l^2 + 14*l + 38/3,
//
// where l is the conic parameter and X=m^2 is already a square.
// A full second-stage point must also satisfy the B and A component
// square conditions.  This script reduces the rank-one Mordell-Weil
// group of E_m modulo good primes and keeps only those residue classes
// whose associated finite-field w passes the A/B tests.
//
// By default BranchMode := "conservative": if the finite reduction hits
// a bad denominator or zero component discriminant, the class is kept at
// that prime.  BranchMode := "strict" rejects zero discriminants; this
// is useful for sieving the nondegenerate open part and for confirming
// that the only obvious survivors in conservative mode are the known
// A-degenerate w = +/-8 branches.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned MaxPrime then MaxPrime := 199; end if;
if Type(MaxPrime) eq MonStgElt then MaxPrime := StringToInteger(MaxPrime); end if;
if not assigned BranchMode then BranchMode := "conservative"; end if;
if not assigned MaxCRTClasses then MaxCRTClasses := 200000; end if;
if Type(MaxCRTClasses) eq MonStgElt then MaxCRTClasses := StringToInteger(MaxCRTClasses); end if;
if not assigned ReportLimit then ReportLimit := 40; end if;
if Type(ReportLimit) eq MonStgElt then ReportLimit := StringToInteger(ReportLimit); end if;

P<lam> := PolynomialRing(Q);
fX := -14*lam^3 - (Q!38/3)*lam^2 + 14*lam + Q!38/3;
C := HyperellipticCurve(fX);
ptsInf := PointsAtInfinity(C);
error if #ptsInf eq 0, "E_m quotient has no rational infinity";
Eraw, phi := EllipticCurve(C, ptsInf[1]);
Emin, mp := MinimalModel(Eraw);
MW, toE := MordellWeilGroup(Emin);
invs := Invariants(MW);
freeidx := [i : i in [1..#invs] | invs[i] eq 0];
torsidx := [i : i in [1..#invs] | invs[i] ne 0];
error if #freeidx ne 1, "expected rank one quotient";

phiInv := Inverse(phi);
mpInv := Inverse(mp);
gensMW := [MW.i : i in [1..#invs]];
freeElt := gensMW[freeidx[1]];
torsElts := [gensMW[i] : i in torsidx];

function RedRat(q, F)
    q := Q!q;
    p := Z!Characteristic(F);
    den := Denominator(q);
    if den mod p eq 0 then
        return false, F!0;
    end if;
    return true, F!Numerator(q)/F!den;
end function;

function PolyReduce(poly, Lp)
    F := BaseRing(Lp);
    out := Lp!0;
    for i in [0..Degree(poly)] do
        ok, c := RedRat(Coefficient(poly, i), F);
        if not ok then
            return false, Lp!0;
        end if;
        out +:= c*(Lp.1)^i;
    end for;
    return true, out;
end function;

function RedCPoint(Cpt, Cfp)
    F := BaseRing(Ambient(Cfp));
    if Cpt[3] eq 0 then
        pts := PointsAtInfinity(Cfp);
        if #pts eq 0 then
            return false, Cfp!0;
        end if;
        return true, pts[1];
    end if;

    okx, xr := RedRat(Cpt[1], F);
    oky, yr := RedRat(Cpt[2], F);
    okz, zr := RedRat(Cpt[3], F);
    if not (okx and oky and okz) then
        return false, Cfp!0;
    end if;
    if zr eq 0 then
        pts := PointsAtInfinity(Cfp);
        if #pts eq 0 then
            return false, Cfp!0;
        end if;
        return true, pts[1];
    end if;
    try
        return true, Cfp![xr, yr, zr];
    catch e
        return false, Cfp!0;
    end try;
end function;

function CPointFromMW(g)
    Pmin := toE(g);
    if Pmin eq Emin!0 then
        return ptsInf[1];
    end if;
    return phiInv(mpInv(Pmin));
end function;

function ReduceMWToFiniteE(g, Cfp, phip)
    Cpt := CPointFromMW(g);
    ok, Cp := RedCPoint(Cpt, Cfp);
    if not ok then
        return false, Codomain(phip)!0;
    end if;
    return true, phip(Cp);
end function;

function TorsionCombos()
    combos := [<MW!0, []>];
    for j in [1..#torsElts] do
        next := [];
        ord := invs[torsidx[j]];
        for c in combos do
            for a in [0..ord-1] do
                Append(~next, <c[1] + a*torsElts[j], c[2] cat [a]>);
            end for;
        end for;
        combos := next;
    end for;
    return combos;
end function;

function CRT2(a, M, b, N)
    a := Z!a; M := Z!M; b := Z!b; N := Z!N;
    g := GCD(M, N);
    delta := b - a;
    if delta mod g ne 0 then
        return false, 0, 0;
    end if;
    l := LCM(M, N);
    M1 := M div g;
    N1 := N div g;
    if N1 eq 1 then
        t := Z!0;
    else
        t := ((delta div g)*InverseMod(M1 mod N1, N1)) mod N1;
    end if;
    return true, (a + M*t) mod l, l;
end function;

function KeyClass(ti, a, M)
    return Sprintf("%o:%o:%o", ti, a, M);
end function;

function SquareNonzero(a)
    return a ne 0 and IsSquare(a);
end function;

function ComponentPassFinite(gq, XR)
    F := BaseRing(Parent(gq));
    if Degree(gq) ne 2 or LeadingCoefficient(gq) eq F!0 then
        return true;  // not a good reduction for this test
    end if;

    dsc := Discriminant(gq);
    if dsc eq F!0 then
        return BranchMode eq "conservative";
    end if;

    if IsSquare(dsc) then
        roots := Roots(gq);
        if #roots lt 2 then
            return true;
        end if;
        for rt in roots do
            val := XR - rt[1];
            if not SquareNonzero(val) then
                return false;
            end if;
        end for;
        return true;
    end if;

    monic := gq/LeadingCoefficient(gq);
    K<th> := ext<F | monic>;
    val := K!XR - th;
    return val ne 0 and IsSquare(val);
end function;

function FullStagePassFinite(wv)
    F := Parent(wv);
    if Characteristic(F) eq 2 then
        return true;
    end if;
    if wv^2 eq F!1 then
        return true;
    end if;

    PF<x> := PolynomialRing(F);
    R := F!-8;
    tv := (2*R^2 + (1-wv^2)*R - 2*wv^2)/(4*(wv^2-1));
    c4v := R + 2 + 4*tv;
    if c4v eq F!0 then
        return true;
    end if;

    Av := x^2 + (R^3 + 4*R^2*tv + R - 8*R*tv + 4*tv)*x + R^4;
    Bv := c4v*x^2 + (R^2 + 4*R + 1 + 8*tv)*x + (2*R^2 + R + 4*tv);
    XR := -c4v*R;
    Atq := PF![co : co in Coefficients(c4v^2*Evaluate(Av, x/c4v))];
    Btq := PF![co : co in Coefficients(c4v*Evaluate(Bv, x/c4v))];

    return ComponentPassFinite(Atq, XR) and ComponentPassFinite(Btq, XR);
end function;

function LambdaToWFinite(lv)
    F := Parent(lv);
    den := lv^2 - F!1;
    if den eq F!0 then
        return false, F!0, F!0;
    end if;
    Xv := (F!-2016*lv - F!1824)/den;
    if Xv eq F!1008 then
        return false, F!0, Xv;
    end if;
    return true, (Xv + F!1008)/(Xv - F!1008), Xv;
end function;

function PrimeAllowedResidues(p, torsCombos)
    F := GF(p);
    Lp<lp> := PolynomialRing(F);
    okf, fXp := PolyReduce(fX, Lp);
    if not okf or Discriminant(fXp) eq F!0 then
        return false, 0, [], 0, 0, 0;
    end if;

    Cfp := HyperellipticCurve(fXp);
    ptsFp := PointsAtInfinity(Cfp);
    if #ptsFp eq 0 then
        return false, 0, [], 0, 0, 0;
    end if;
    Efp, phip := EllipticCurve(Cfp, ptsFp[1]);

    okG, Gp := ReduceMWToFiniteE(freeElt, Cfp, phip);
    if not okG then
        return false, 0, [], 0, 0, 0;
    end if;
    ordG := Order(Gp);
    if ordG eq 0 then
        return false, 0, [], 0, 0, 0;
    end if;

    redTors := [];
    for tc in torsCombos do
        okT, Tp := ReduceMWToFiniteE(tc[1], Cfp, phip);
        if not okT then
            return false, 0, [], 0, 0, 0;
        end if;
        Append(~redTors, Tp);
    end for;

    lambdaByE := AssociativeArray();
    finiteLambdaCount := 0;
    for Cp in Points(Cfp) do
        Ep := phip(Cp);
        if Cp[3] ne 0 then
            lambdaByE[Ep] := Cp[1]/Cp[3];
            finiteLambdaCount +:= 1;
        end if;
    end for;

    allowedByT := [];
    allowedTotal := 0;
    for ti in [1..#torsCombos] do
        allowed := {};
        for k in [0..ordG-1] do
            Qp := redTors[ti] + k*Gp;
            if not IsDefined(lambdaByE, Qp) then
                Include(~allowed, k);
                continue;
            end if;
            lv := lambdaByE[Qp];
            okw, wv, Xv := LambdaToWFinite(lv);
            if not okw then
                Include(~allowed, k);
                continue;
            end if;
            if FullStagePassFinite(wv) then
                Include(~allowed, k);
            end if;
        end for;
        allowedTotal +:= #allowed;
        Append(~allowedByT, allowed);
    end for;

    return true, ordG, allowedByT, #Efp, finiteLambdaCount, allowedTotal;
end function;

torsCombos := TorsionCombos();
classes := [<ti, 0, 1> : ti in [1..#torsCombos]];
usedPrimes := [];
truncated := false;

printf "R=-8 full-stage MW sieve on E_m quotient\n";
printf "BranchMode=%o MaxPrime=%o MaxCRTClasses=%o\n", BranchMode, MaxPrime, MaxCRTClasses;
printf "fX=%o\n", fX;
printf "Eraw=%o\n", Eraw;
printf "Emin=%o\n", Emin;
printf "RankBounds=%o MW invariants=%o torsion combos=%o\n", RankBounds(Emin), invs, #torsCombos;
printf "initial classes=%o\n", #classes;

for p in PrimesUpTo(MaxPrime) do
    if p in {2,3,7} then
        continue;
    end if;

    ok, ordG, allowedByT, ecount, lambdaCount, allowedTotal := PrimeAllowedResidues(p, torsCombos);
    if not ok then
        printf "p=%o skipped\n", p;
        continue;
    end if;
    Append(~usedPrimes, p);

    before := #classes;
    next := [];
    seen := {};
    for cls in classes do
        ti := cls[1];
        for k in allowedByT[ti] do
            okcrt, r, M := CRT2(cls[2], cls[3], k, ordG);
            if okcrt then
                key := KeyClass(ti, r, M);
                if key notin seen then
                    Include(~seen, key);
                    Append(~next, <ti, r, M>);
                    if #next gt MaxCRTClasses then
                        truncated := true;
                        break cls;
                    end if;
                end if;
            end if;
        end for;
    end for;
    classes := next;

    printf "p=%o #E(Fp)=%o ordG=%o finite_lambdas=%o allowed_residue_sum=%o classes %o -> %o",
        p, ecount, ordG, lambdaCount, allowedTotal, before, #classes;
    if truncated then
        printf " TRUNCATED";
    end if;
    printf "\n";
    if #classes le ReportLimit and not truncated then
        printf "  survivors=%o\n", classes;
    end if;
    if #classes eq 0 or truncated then
        break;
    end if;
end for;

printf "DONE\n";
printf "used_primes=%o\n", usedPrimes;
printf "final_survivors=%o truncated=%o\n", #classes, truncated;
if #classes le ReportLimit and not truncated then
    printf "final_classes=%o\n", classes;
end if;

quit;
