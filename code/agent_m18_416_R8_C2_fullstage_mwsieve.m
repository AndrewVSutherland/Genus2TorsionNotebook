//////////////////////////////////////////////////////////////////////
// Mordell-Weil sieve on the original R = -8 C2 elliptic fiber.
//
// This is sharper than sieving only on the S_B lambda quotient: a point
// of the C2 fiber carries both m and the chosen square-root of
// G(R,m).  We reduce the rank-one Mordell-Weil group of this elliptic
// curve modulo good primes and retain only residue classes for which
// the finite-field reductions satisfy both second-stage square tests
// S_A and S_B.
//
// BranchMode := "conservative" keeps classes at bad denominators or
// zero reduced square-values.  BranchMode := "strict" rejects zero
// square-values, i.e. it sieves the nondegenerate open part.
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

// ---- R = -8 symbolic second-stage functions over Q(m) ----
K1<m> := RationalFunctionField(Q);
PX<x> := PolynomialRing(K1);
RvK := K1!-8;
KvK := -2*RvK*(RvK^2-1);
wK := (m^2 + KvK)/(m^2 - KvK);
tK := (2*RvK^2 + (1-wK^2)*RvK - 2*wK^2)/(4*(wK^2-1));
c4K := RvK + 2 + 4*tK;
Apol := x^2 + (RvK^3 + 4*RvK^2*tK + RvK - 8*RvK*tK + 4*tK)*x + RvK^4;
Bpol := c4K*x^2 + (RvK^2 + 4*RvK + 1 + 8*tK)*x + (2*RvK^2 + RvK + 4*tK);
XRK := -c4K*RvK;
At := PX![K1!co : co in Coefficients(c4K^2*Evaluate(Apol, x/c4K))];
Bt := PX![K1!co : co in Coefficients(c4K*Evaluate(Bpol, x/c4K))];
GK := 2*(RvK^2-1)*(RvK*(2*RvK+1) - wK^2*(RvK+2));
alphaA := XRK + Coefficient(At,1)/2;
alphaB := XRK + Coefficient(Bt,1)/2;
okA, hA := IsSquare(Evaluate(At, XRK)/GK); assert okA;
okB, hB := IsSquare(Evaluate(Bt, XRK)/GK); assert okB;

P<z> := PolynomialRing(Q);
Rv := Q!-8;
Kv := -2*Rv*(Rv^2-1);
Gam := 2*(Rv^2-1)*(Rv*(2*Rv+1)*(z^2-Kv)^2 - (Rv+2)*(z^2+Kv)^2);
C := HyperellipticCurve(Gam);
m0 := Q!-28;
okg0, g0 := IsSquare(Evaluate(Gam, m0));
error if not okg0, "base point is not on C2 fiber";
P0 := C![m0, g0, Q!1];
Eraw, phi := EllipticCurve(C, P0);
Emin, mp := MinimalModel(Eraw);
MW, toE := MordellWeilGroup(Emin);
invs := Invariants(MW);
freeidx := [i : i in [1..#invs] | invs[i] eq 0];
torsidx := [i : i in [1..#invs] | invs[i] ne 0];
error if #freeidx ne 1, "expected rank one C2 fiber";

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

function EvalPolyQAt(poly, val)
    F := Parent(val);
    out := F!0;
    for i in [0..Degree(poly)] do
        ok, c := RedRat(Coefficient(poly, i), F);
        if not ok then
            return false, F!0;
        end if;
        out +:= c*val^i;
    end for;
    return true, out;
end function;

function EvalRatK1At(rat, val)
    F := Parent(val);
    num := Numerator(rat);
    den := Denominator(rat);
    okn, nv := EvalPolyQAt(num, val);
    okd, dv := EvalPolyQAt(den, val);
    if not (okn and okd) or dv eq F!0 then
        return false, F!0;
    end if;
    return true, nv/dv;
end function;

function RedCPoint(Cpt, Cfp)
    F := BaseRing(Ambient(Cfp));
    okx, xr := RedRat(Cpt[1], F);
    oky, yr := RedRat(Cpt[2], F);
    okz, zr := RedRat(Cpt[3], F);
    if not (okx and oky and okz) then
        return false, Cfp!0;
    end if;
    try
        return true, Cfp![xr, yr, zr];
    catch e
        return false, Cfp!0;
    end try;
end function;

function CPointFromMW(g)
    return phiInv(mpInv(toE(g)));
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

function SquareOK(a)
    if BranchMode eq "conservative" then
        return a eq 0 or IsSquare(a);
    end if;
    return a ne 0 and IsSquare(a);
end function;

function FullStagePassAtC2Point(mv, gv)
    F := Parent(mv);
    if Characteristic(F) in {2,3,7} then
        return true;
    end if;
    denG := mv^2 - F!1008;
    if denG eq F!0 then
        return true;
    end if;
    sg := gv/denG; // square-root of G(R,m), since gv^2 = G*(m^2-K)^2

    oka, aA := EvalRatK1At(alphaA, mv);
    okh, hAv := EvalRatK1At(hA, mv);
    okb, aB := EvalRatK1At(alphaB, mv);
    okhb, hBv := EvalRatK1At(hB, mv);
    if not (oka and okh and okb and okhb) then
        return true;
    end if;

    passA := SquareOK(2*aA + 2*hAv*sg) or SquareOK(2*aA - 2*hAv*sg);
    if not passA then
        return false;
    end if;
    passB := SquareOK(2*aB + 2*hBv*sg) or SquareOK(2*aB - 2*hBv*sg);
    return passB;
end function;

function PrimeAllowedResidues(p, torsCombos)
    F := GF(p);
    Lp<zp> := PolynomialRing(F);
    okf, Gamp := PolyReduce(Gam, Lp);
    if not okf or Discriminant(Gamp) eq F!0 then
        return false, 0, [], 0, 0, 0;
    end if;

    Cfp := HyperellipticCurve(Gamp);
    okP0, P0p := RedCPoint(P0, Cfp);
    if not okP0 then
        return false, 0, [], 0, 0, 0;
    end if;
    Efp, phip := EllipticCurve(Cfp, P0p);

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

    cByE := AssociativeArray();
    finiteCount := 0;
    for Cp in Points(Cfp) do
        Ep := phip(Cp);
        if Cp[3] ne F!0 then
            cByE[Ep] := <Cp[1]/Cp[3], Cp[2]/Cp[3]^2>;
            finiteCount +:= 1;
        end if;
    end for;

    allowedByT := [];
    allowedTotal := 0;
    for ti in [1..#torsCombos] do
        allowed := {};
        for k in [0..ordG-1] do
            Qp := redTors[ti] + k*Gp;
            if not IsDefined(cByE, Qp) then
                Include(~allowed, k);
                continue;
            end if;
            cg := cByE[Qp];
            if FullStagePassAtC2Point(cg[1], cg[2]) then
                Include(~allowed, k);
            end if;
        end for;
        allowedTotal +:= #allowed;
        Append(~allowedByT, allowed);
    end for;

    return true, ordG, allowedByT, #Efp, finiteCount, allowedTotal;
end function;

torsCombos := TorsionCombos();
classes := [<ti, 0, 1> : ti in [1..#torsCombos]];
usedPrimes := [];
truncated := false;

printf "R=-8 C2 full-stage MW sieve\n";
printf "BranchMode=%o MaxPrime=%o MaxCRTClasses=%o\n", BranchMode, MaxPrime, MaxCRTClasses;
printf "Gam=%o\n", Gam;
printf "Eraw=%o\n", Eraw;
printf "Emin=%o\n", Emin;
printf "RankBounds=%o MW invariants=%o torsion combos=%o\n", RankBounds(Emin), invs, #torsCombos;
printf "initial classes=%o\n", #classes;

for p in PrimesUpTo(MaxPrime) do
    if p in {2,3,7} then
        continue;
    end if;

    ok, ordG, allowedByT, ecount, finiteCount, allowedTotal := PrimeAllowedResidues(p, torsCombos);
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

    printf "p=%o #E(Fp)=%o ordG=%o finite_C_points=%o allowed_residue_sum=%o classes %o -> %o",
        p, ecount, ordG, finiteCount, allowedTotal, before, #classes;
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
