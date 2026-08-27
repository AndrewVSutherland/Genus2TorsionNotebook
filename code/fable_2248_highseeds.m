//////////////////////////////////////////////////////////////////////
// fable_2248_highseeds.m  (2026-07-18, Fable session)
//
// Depth-1 Richelot fans of the SEVEN highest-value [2,2,2,8] seeds:
// the six certified-simple off-rectangle curves at heights 1e13-1e16
// (already beyond the d<=65535 audit box) and the twisted
// (29,121,125,145) curve.  Model: y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2).
// SRCCHK verifies the model convention via 64 | #J(F_p) on the source.
//
// Run:  magma -b code/fable_2248_highseeds.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    den := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(den^2*f), den;
end function;

function NormalizePoly(D)
    fD, hD := HyperellipticPolynomials(D);
    FD := hD eq 0 select P!fD else P!(hD^2 + 4*fD);
    FI, _ := IntegralSquareScale(FD);
    return FI;
end function;

function GoodOrders(f, nmax)   // up to nmax good-prime #J(F_p) values
    Zf := P![Z!c : c in Coefficients(f)];
    orders := [];
    for p in [13,17,19,23,29,31,37,41,43,47] do
        if #orders ge nmax then break; end if;
        Fp := GF(p); Pp := PolynomialRing(Fp);
        fp := Pp![Fp!c : c in Coefficients(Zf)];
        if Degree(fp) ne Degree(Zf) then continue; end if;
        if not IsSeparable(fp) then continue; end if;
        okp := true;
        try
            Np := Evaluate(LPolynomial(HyperellipticCurve(fp)), 1);
        catch ee okp := false; end try;
        if okp then Append(~orders, Np); end if;
    end for;
    return orders;
end function;

function PassGate(f)
    os := GoodOrders(f, 4);
    if #os lt 2 then return false; end if;
    g := 0;
    for Np in os do g := Gcd(g, Np); end for;
    return g mod 96 eq 0;
end function;

procedure FanSource(f, tag, ~nGate, ~nBig)
    ok := true;
    J := 0;
    try
        J := Jacobian(HyperellipticCurve(f));
    catch ee ok := false; end try;
    if not ok then printf "BADSOURCE %o\n", tag; return; end if;
    Rs := [];
    try
        Rs := RichelotIsogenousSurfaces(J);
    catch ee
        printf "RICHELOT_FAIL %o\n", tag; return;
    end try;
    printf "FAN %o codomains=%o\n", tag, #Rs;
    for i in [1..#Rs] do
        S := Rs[i];
        tS := Type(S);
        if not (tS eq JacHyp or tS eq CrvHyp) then
            printf "NONJAC %o.%o type=%o\n", tag, i, tS;
            continue;
        end if;
        D := tS eq JacHyp select Curve(S) else S;
        fN := NormalizePoly(D);
        if not PassGate(fN) then continue; end if;
        nGate +:= 1;
        okT := true;
        try
            T := TorsionSubgroup(Jacobian(HyperellipticCurve(fN)));
        catch ee okT := false; end try;
        if not okT then printf "TORSFAIL %o.%o\n", tag, i; continue; end if;
        inv := Invariants(T);
        ord := #inv eq 0 select 1 else &*inv;
        if inv eq [2,2,4,8] then
            printf "JACKPOT %o.%o inv=%o f=%o\n", tag, i, inv, fN;
        elif ord ge 96 then
            nBig +:= 1;
            printf "BIG %o.%o inv=%o order=%o f=%o\n", tag, i, inv, ord, fN;
        else
            printf "GATEPASS %o.%o inv=%o\n", tag, i, inv;
        end if;
    end for;
end procedure;

tuples := [
    [69780688487820, 1325833081268580, 2946719342448225, 3981345244907824],
    [87225860609775, 117851829446096, 261930608217620, 4976681556134780],
    [3795310132322, 4664593520206, 29610028432612, 3186162856084319],
    [405616827844, 43646066521703, 277057639659506, 340515326975038],
    [42539545655911, 150582462498800, 662625094557200, 74876635684963600],
    [376456156247, 42539545655911, 187191589212409, 662625094557200],
    [29, 121, 125, 145]
];

ng := 0; nb := 0;
for k in [1..#tuples] do
    t := tuples[k];
    f := x*(x+t[1]^2)*(x+t[2]^2)*(x+t[3]^2)*(x+t[4]^2);
    os := GoodOrders(f, 2);
    srcok := #os ge 1 and &and[Np mod 64 eq 0 : Np in os];
    printf "SRCCHK H%o model64=%o\n", k, srcok;
    if not srcok then
        printf "SKIP H%o (model convention mismatch?)\n", k;
        continue;
    end if;
    FanSource(f, Sprintf("H%o", k), ~ng, ~nb);
end for;
printf "HIGHSEEDS_DONE gate=%o big=%o\n", ng, nb;
quit;
