//////////////////////////////////////////////////////////////////////
// fable_2248_bank_walk.m  (2026-07-18, Fable session)
//
// Wide shallow Richelot walk hunting [2,2,4,8] (order 128), and any
// torsion order >= 96, among depth-1 Richelot codomains of the
// project's (2,2,4,4) and (2,2,2,8) tuple banks.
//
// Why: direct all-squares search is closed for d<=65535 (prod-07
// sign-reduction + tier-1 audit); but codomain torsion is NOT pulled
// back through the isogeny (the record's own [2,12] leaves sit one
// Richelot step from order 96), so walking from bank members samples
// high-height curves the audit cannot reach.  Simplicity is isogeny-
// invariant, so codomains of certified banks need no new certificate
// for the *search* phase (hits get the full checklist separately).
//
// Funnel per codomain: normalize -> two good primes from {13,...,37}
// with 96 | #J(F_p) -> exact TorsionSubgroup.  [2,2,4,8] would give
// 128 | #J(F_p); the 96 gate also catches any new order>=96 group.
//
// Stage 0 validates the funnel end-to-end on the record's graph.
// Stage 1 control: depth-1 fan of the known SPLIT [2,2,4,8] tuple.
// Stage 2 production: stride samples of tor2244.txt and tor2228.txt.
//
// Run (from repo root):
//   magma -b N44:=150 N28:=150 code/fable_2248_bank_walk.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

if not assigned N44 then N44 := 150;
elif Type(N44) eq MonStgElt then N44 := StringToInteger(N44); end if;
if not assigned N28 then N28 := 150;
elif Type(N28) eq MonStgElt then N28 := StringToInteger(N28); end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

DATA := "paper/scripts_and_data/";

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

// gcd order gate over up to 4 good primes: 96 | gcd(#J(F_p)) required.
function PassGate(f)
    Zf := P![Z!c : c in Coefficients(f)];
    g := 0; good := 0;
    for p in [13,17,19,23,29,31,37,41] do
        if good ge 4 then break; end if;
        Fp := GF(p);
        Pp := PolynomialRing(Fp);
        fp := Pp![Fp!c : c in Coefficients(Zf)];
        if Degree(fp) ne Degree(Zf) then continue; end if;
        if not IsSeparable(fp) then continue; end if;
        okp := true;
        try
            Cp := HyperellipticCurve(fp);
            Np := Evaluate(LPolynomial(Cp), 1);
        catch ee okp := false; end try;
        if not okp then continue; end if;
        good +:= 1;
        g := Gcd(g, Np);
        if g mod 96 ne 0 then return false; end if;   // early kill
    end for;
    return good ge 2 and g mod 96 eq 0;
end function;

// full fan of one source polynomial
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
    for i in [1..#Rs] do
        S := Rs[i];
        tS := Type(S);
        if not (tS eq JacHyp or tS eq CrvHyp) then continue; end if;
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

// ---------- Stage 0: validation on the record graph ----------
frec := 3027600*x^6 + 2950382280*x^5 + 602288814361*x^4
        - 63417934304484*x^3 - 2122595910966478*x^2
        + 128056619498204124*x + 3322970988364151397;
Jr := Jacobian(HyperellipticCurve(frec));
Rr := RichelotIsogenousSurfaces(Jr);
leaf := 0;
for S in Rr do
    tS := Type(S);
    if tS eq JacHyp or tS eq CrvHyp then
        leaf := NormalizePoly(tS eq JacHyp select Curve(S) else S);
        break;
    end if;
end for;
ng := 0; nb := 0;
FanSource(leaf, "VALIDATE", ~ng, ~nb);
if nb ge 1 then
    printf "VALIDATE_OK big_hits=%o (record rediscovered from its leaf)\n", nb;
else
    printf "VALIDATE_FAIL — funnel did not rediscover the record; DO NOT TRUST THIS RUN\n";
end if;

// ---------- Stage 1: control (known split [2,2,4,8] tuple) ----------
lines := Split(Read(DATA cat "ten2248models_abcd.txt"), "\n");
ctrl := eval lines[1];
a := ctrl[1]; b := ctrl[2]; c := ctrl[3]; d := ctrl[4];
fctl := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
ng0 := 0; nb0 := 0;
FanSource(fctl, "CONTROL", ~ng0, ~nb0);
printf "CONTROL_DONE gatepass=%o big=%o\n", ng0, nb0;

// ---------- Stage 2: bank walks ----------
procedure WalkBank(fname, N, tagpre, ~tGate, ~tBig)
    lines := Split(Read(DATA cat fname), "\n");
    n := #lines;
    step := Maximum(1, n div N);
    done := 0;
    i := 1;
    while i le n and done lt N do
        line := lines[i];
        if #line gt 2 then
            tup := eval line;
            aa := tup[1]; bb := tup[2]; cc := tup[3]; dd := tup[4];
            f := x*(x+aa^2)*(x+bb^2)*(x+cc^2)*(x+dd^2);
            FanSource(f, Sprintf("%o_r%o", tagpre, i), ~tGate, ~tBig);
            done +:= 1;
            if done mod 10 eq 0 then
                printf "PROGRESS %o done=%o row=%o gate=%o big=%o\n", tagpre, done, i, tGate, tBig;
            end if;
        end if;
        i +:= step;
    end while;
    printf "BANK_DONE %o sources=%o\n", tagpre, done;
end procedure;

g44 := 0; b44 := 0;
WalkBank("tor2244.txt", N44, "B44", ~g44, ~b44);
g28 := 0; b28 := 0;
WalkBank("tor2228.txt", N28, "B28", ~g28, ~b28);

printf "SEARCH_DONE gate44=%o big44=%o gate28=%o big28=%o\n", g44, b44, g28, b28;
quit;
