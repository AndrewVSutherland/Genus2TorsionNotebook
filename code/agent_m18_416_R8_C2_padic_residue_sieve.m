//////////////////////////////////////////////////////////////////////
// p-adic residue version of the R = -8 C2 Mordell-Weil sieve.
//
// For a chosen prime p, enumerate the rank-one MW classes
//     T + n*G,  n modulo ord(red(G))*p^Depth,
// map them back to the C2 quartic, and test the actual Q_p second-stage
// S_A and S_B square conditions.  This is not a formal power-series
// Chabauty proof by itself, but it identifies whether higher p-adic
// residue information can plausibly close the remaining rank-one fiber.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned Prime then Prime := 11; end if;
if Type(Prime) eq MonStgElt then Prime := StringToInteger(Prime); end if;
if not assigned Depth then Depth := 2; end if;
if Type(Depth) eq MonStgElt then Depth := StringToInteger(Depth); end if;
if not assigned Prec then Prec := 80; end if;
if Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
if not assigned ReportLimit then ReportLimit := 30; end if;
if Type(ReportLimit) eq MonStgElt then ReportLimit := StringToInteger(ReportLimit); end if;

K2<R,m> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);
K := -2*R*(R^2-1);
w := (m^2 + K)/(m^2 - K);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
c4 := R + 2 + 4*t;
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := c4*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
XR := -c4*R;
At := PX![K2!co : co in Coefficients(c4^2*Evaluate(Apol, x/c4))];
Bt := PX![K2!co : co in Coefficients(c4*Evaluate(Bpol, x/c4))];
Gfun := 2*(R^2-1)*(R*(2*R+1) - w^2*(R+2));
alphaA := XR + Coefficient(At,1)/2;
alphaB := XR + Coefficient(Bt,1)/2;
okA, hA := IsSquare(Evaluate(At, XR)/Gfun); assert okA;
okB, hB := IsSquare(Evaluate(Bt, XR)/Gfun); assert okB;
funcs := [Gfun, alphaA, hA, alphaB, hB];
fnum := [Numerator(f) : f in funcs];
fden := [Denominator(f) : f in funcs];

Rv := Q!-8;
P<z> := PolynomialRing(Q);
Kv := -2*Rv*(Rv^2-1);
Gam := 2*(Rv^2-1)*(Rv*(2*Rv+1)*(z^2-Kv)^2 - (Rv+2)*(z^2+Kv)^2);
C := HyperellipticCurve(Gam);
m0 := Q!-28;
okg0, g0 := IsSquare(Evaluate(Gam, m0));
assert okg0;
P0 := C![m0, g0, Q!1];
Eraw, phi := EllipticCurve(C, P0);
Emin, mp := MinimalModel(Eraw);
MW, toE := MordellWeilGroup(Emin);
invs := Invariants(MW);
freeidx := [i : i in [1..#invs] | invs[i] eq 0];
torsidx := [i : i in [1..#invs] | invs[i] ne 0];
assert #freeidx eq 1;
gensMW := [MW.i : i in [1..#invs]];
Gmw := gensMW[freeidx[1]];
torsElts := [gensMW[i] : i in torsidx];
phiInv := Inverse(phi);
mpInv := Inverse(mp);

function EvalAll(mv)
    vals := [];
    for i in [1..#funcs] do
        dd := Evaluate(fden[i], [Rv, mv]);
        if dd eq 0 or Valuation(dd) ne 0 then
            return false, [];
        end if;
        Append(~vals, Evaluate(fnum[i], [Rv, mv])/dd);
    end for;
    return true, vals;
end function;

function PointPasses(mv, Qp)
    ok, v := EvalAll(mv);
    if not ok then
        return true;
    end if;
    Gv := v[1]; aA := v[2]; hAv := v[3]; aB := v[4]; hBv := v[5];
    if Gv eq 0 or hAv eq 0 or hBv eq 0 then
        return true;
    end if;
    okG, sg := IsSquare(Qp!Gv);
    if not okG then
        return false;
    end if;
    passA := false;
    for sgn in [1,-1] do
        VA := Qp!(2*aA) + sgn*(Qp!(2*hAv))*sg;
        if VA ne 0 and IsSquare(VA) then
            passA := true;
            break;
        end if;
    end for;
    if not passA then
        return false;
    end if;
    for sgn in [1,-1] do
        VB := Qp!(2*aB) + sgn*(Qp!(2*hBv))*sg;
        if VB ne 0 and IsSquare(VB) then
            return true;
        end if;
    end for;
    return false;
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

function RedRat(q, F)
    q := Q!q;
    den := Denominator(q);
    if den mod Characteristic(F) eq 0 then
        return false, F!0;
    end if;
    return true, F!Numerator(q)/F!den;
end function;

F := GF(Prime);
Ered := EllipticCurve([F!a : a in aInvariants(Emin)]);
Gpt := toE(Gmw);
if Gpt eq Emin!0 then
    error "free generator maps to zero";
end if;
okx, gx := RedRat(Gpt[1]/Gpt[3], F);
oky, gy := RedRat(Gpt[2]/Gpt[3], F);
error if not (okx and oky), "generator has bad reduction";
Gred := Ered![gx, gy, F!1];
ordRed := Order(Gred);
modulus := ordRed * Prime^Depth;
Qp := pAdicField(Prime, Prec);
Eminp := BaseChange(Emin, Qp);
GptMinQ := toE(Gmw);
GptMin := Eminp![Qp!(GptMinQ[1]/GptMinQ[3]), Qp!(GptMinQ[2]/GptMinQ[3]), Qp!1];
torsCombos := TorsionCombos();

printf "R=-8 p-adic MW residue sieve, p=%o Depth=%o Prec=%o\n", Prime, Depth, Prec;
printf "Emin=%o MW invariants=%o ordRed(G)=%o modulus=%o torsion_combos=%o\n",
    Emin, invs, ordRed, modulus, #torsCombos;

survivors := [];
tested := 0;
badmap := 0;
for ti in [1..#torsCombos] do
    Telt := torsCombos[ti][1];
    TptQ := toE(Telt);
    if TptQ eq Emin!0 then
        Pmin := Eminp!0;
    else
        Pmin := Eminp![Qp!(TptQ[1]/TptQ[3]), Qp!(TptQ[2]/TptQ[3]), Qp!1];
    end if;
    for n in [0..modulus-1] do
        tested +:= 1;
        if Pmin eq Eminp!0 then
            mv := Qp!m0;
        else
            xx := Pmin[1]/Pmin[3];
            yy := Pmin[2]/Pmin[3];
            denm := 2*xx + yy - 128;
            if denm eq 0 or Valuation(denm) ne 0 then
                badmap +:= 1;
                Append(~survivors, <ti,n,"mapden">);
                Pmin +:= GptMin;
                continue;
            end if;
            try
                mv := -28*(6*xx + yy - 24)/denm;
            catch e
                badmap +:= 1;
                Append(~survivors, <ti,n,"mapden">);
                Pmin +:= GptMin;
                continue;
            end try;
        end if;
        if PointPasses(mv, Qp) then
            Append(~survivors, <ti,n,"pass">);
        end if;
/*
        catch e
            badmap +:= 1;
            Append(~survivors, <ti,n,"mapfail">);
        end try;
*/
        Pmin +:= GptMin;
    end for;
    printf "after torsion combo %o survivors=%o\n", ti, #survivors;
end for;

printf "tested=%o badmap_or_infinity=%o survivors=%o\n", tested, badmap, #survivors;
if #survivors le ReportLimit then
    printf "survivors=%o\n", survivors;
else
    printf "first_survivors=%o\n", survivors[1..ReportLimit];
end if;
print "DONE";
quit;
