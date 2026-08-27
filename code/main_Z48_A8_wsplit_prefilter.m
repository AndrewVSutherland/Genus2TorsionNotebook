//////////////////////////////////////////////////////////////////////
//  Z/48 scout on the W-split sublocus of the A(8) chart.
//
//  This adapts the W-split [4,16]/[2,2,16] search to the Z/48 question:
//  first require 48 | #J(F_p) at several good primes, then test exact
//  halving of the visible order-8 class.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if assigned RH and Type(RH) eq MonStgElt then RH := StringToInteger(RH); end if;
if not assigned RH then RH := 6; end if;
if assigned TH and Type(TH) eq MonStgElt then TH := StringToInteger(TH); end if;
if not assigned TH then TH := 6; end if;
if assigned BH and Type(BH) eq MonStgElt then BH := StringToInteger(BH); end if;
if not assigned BH then BH := 6; end if;
if assigned Progress and Type(Progress) eq MonStgElt then Progress := StringToInteger(Progress); end if;
if not assigned Progress then Progress := 20000; end if;
if assigned MinGood and Type(MinGood) eq MonStgElt then MinGood := StringToInteger(MinGood); end if;
if not assigned MinGood then MinGood := 2; end if;
if assigned SkipEven and Type(SkipEven) eq MonStgElt then SkipEven := SkipEven in {"true","True","1","yes","Yes"}; end if;
if not assigned SkipEven then SkipEven := true; end if;
if assigned PositiveBetaOnly and Type(PositiveBetaOnly) eq MonStgElt then
    PositiveBetaOnly := PositiveBetaOnly in {"true","True","1","yes","Yes"};
end if;
if not assigned PositiveBetaOnly then PositiveBetaOnly := true; end if;

filter_primes := [5,7,11,13,17,19,23,29,31];

function HeightRationals(H)
    vals := [Q!0];
    for den in [1..H] do
        for num in [-H..H] do
            if GCD(num, den) eq 1 then
                Append(~vals, Q!num/Q!den);
            end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

function IntegralSquareModel(f)
    L := 1;
    for c in Coefficients(f) do
        L := LCM(L, Denominator(c));
    end for;
    return P!(L^2*f);
end function;

function IsEvenPolynomial(f)
    for i in [0..Degree(f)] do
        if IsOdd(i) and Coefficient(f, i) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function A8f(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv;
    d := e + 2*pv - rv^2;
    lambda := rv/tv;
    u := pv + rv*tv - 2*rv;
    v := e + rv^2 - rv*pv - rv^2*tv + 3*pv*tv - rv*tv^2;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d;
    q := a*x^2 + b*x + c;
    L := rv*x + (pv - rv^2);
    g8 := x^2 + u*x + v;
    f := q*(Qpoly^2 + q);
    ellBase := -(q + Qpoly*L);
    return f, g8, ellBase;
end function;

function SplitWPvals(rv, tv, bv)
    Pp<pp> := PolynomialRing(Q);
    e := tv^2 - 2*pp*tv/rv;
    d := e + 2*pp - rv^2;
    lambda := rv/tv;
    a := rv^2 - lambda;
    b := 2*rv*pp - 2*lambda*(pp + rv*tv) + 2*rv*lambda;
    c := pp^2 + 2*pp*rv^2 - rv^4 - rv^3*tv - rv*pp^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pp + rv^2*tv - 3*pp*tv + rv*tv^2);
    Pc := 2*d + a;
    Q1 := b;
    Rc := d^2 + c;
    E := bv^2*(Pc + bv^2)^2 - 4*Rc*bv^2 - Q1^2;
    En := Numerator(E);
    if En eq 0 or Degree(En) lt 1 then
        return [];
    end if;
    return [rt[1] : rt in Roots(En)];
end function;

function SquareQuarticHalves(u, v, f)
    A<M,N> := PolynomialRing(Q, 2);
    RX<X> := PolynomialRing(A);
    phi := hom<P -> RX | X>;
    uX := phi(u);
    fX := phi(f);
    ell1 := phi(v) + uX*(M*X + N);
    if (ell1^2 - fX) mod uX ne 0 then
        return [];
    end if;
    S := ExactQuotient(ell1^2 - fX, uX);
    if Degree(S) ne 4 then
        return [];
    end if;
    s4 := Coefficient(S,4); s3 := Coefficient(S,3); s2 := Coefficient(S,2);
    s1 := Coefficient(S,1); s0 := Coefficient(S,0);
    E1 := 8*s4^2*s1 - s3*(4*s4*s2 - s3^2);
    E0 := 64*s4^3*s0 - (4*s4*s2 - s3^2)^2;
    try
        pts := Variety(ideal<A | E1, E0>);
    catch e
        return [];
    end try;
    out := [];
    for pt in pts do
        Mv := Q!pt[1];
        Nv := Q!pt[2];
        ell1q := v + u*(Mv*x + Nv);
        if (ell1q^2 - f) mod u ne 0 then
            continue;
        end if;
        Sq := ExactQuotient(ell1q^2 - f, u);
        if Degree(Sq) ne 4 then
            continue;
        end if;
        s4q := Coefficient(Sq,4);
        if s4q eq 0 then
            continue;
        end if;
        s3q := Coefficient(Sq,3);
        s2q := Coefficient(Sq,2);
        G := x^2 + (s3q/(2*s4q))*x + (4*s4q*s2q - s3q^2)/(8*s4q^2);
        Append(~out, <G, ell1q>);
    end for;
    return out;
end function;

function GoodModelAt(f, p)
    F := GF(p);
    PF<xp> := PolynomialRing(F);
    fp := PF![F!Coefficient(f, i) : i in [0..Degree(f)]];
    if Degree(fp) lt 5 or Discriminant(fp) eq 0 then
        return false, fp;
    end if;
    return true, fp;
end function;

function Passes48Filter(f)
    good := 0;
    used := [];
    for p in filter_primes do
        ok, fp := GoodModelAt(f, p);
        if not ok then
            continue;
        end if;
        n := Integers()!#Jacobian(HyperellipticCurve(fp));
        Append(~used, <p,n>);
        good +:= 1;
        if n mod 48 ne 0 then
            return false, good, used, p, n;
        end if;
        if good ge MinGood then
            return true, good, used, 0, 0;
        end if;
    end for;
    return false, good, used, 0, 0;
end function;

print "Z48_A8_WSPLIT_PREFILTER";
printf "RH=%o TH=%o BH=%o MinGood=%o SkipEven=%o PositiveBetaOnly=%o filter_primes=%o\n",
    RH, TH, BH, MinGood, SkipEven, PositiveBetaOnly, filter_primes;

rvals := HeightRationals(RH);
tvals := HeightRationals(TH);
bvals := HeightRationals(BH);

tested := 0;
pvalsFound := 0;
smooth := 0;
evenSkip := 0;
factorRich := 0;
filterPass := 0;
d8ok := 0;
halveTried := 0;
halvable := 0;
z48Hits := 0;
firstKills := AssociativeArray();

for rv in rvals do
    if rv eq 0 or rv eq 1 then
        continue;
    end if;
    for tv in tvals do
        if tv eq 0 or rv*tv eq 1 then
            continue;
        end if;
        for bv in bvals do
            if bv eq 0 then
                continue;
            end if;
            if PositiveBetaOnly and bv lt 0 then
                continue;
            end if;
            tested +:= 1;
            if Progress gt 0 and tested mod Progress eq 0 then
                printf "PROGRESS tested=%o pvals=%o smooth=%o evenSkip=%o factorRich=%o filterPass=%o d8ok=%o halvable=%o z48Hits=%o\n",
                    tested, pvalsFound, smooth, evenSkip, factorRich, filterPass, d8ok, halvable, z48Hits;
            end if;

            pvs := [];
            try
                pvs := SplitWPvals(rv, tv, bv);
            catch e
                continue;
            end try;
            for pv in pvs do
                pvalsFound +:= 1;
                f, g8, ellBase := A8f(rv, pv, tv);
                if Degree(f) lt 5 or Discriminant(f) eq 0 then
                    continue;
                end if;
                smooth +:= 1;
                if SkipEven and IsEvenPolynomial(f) then
                    evenSkip +:= 1;
                    continue;
                end if;
                if #Factorization(f) lt 3 then
                    continue;
                end if;
                factorRich +:= 1;
                fInt := IntegralSquareModel(f);
                pass48, goodCount, used, killp, killn := Passes48Filter(fInt);
                if not pass48 then
                    key := (killp eq 0) select "insufficient_good" else IntegerToString(killp);
                    if IsDefined(firstKills, key) then
                        firstKills[key] +:= 1;
                    else
                        firstKills[key] := 1;
                    end if;
                    continue;
                end if;
                filterPass +:= 1;

                C := HyperellipticCurve(f);
                J := Jacobian(C);
                ZJ := J!0;
                v8 := (-ellBase) mod g8;
                D8 := J![g8, v8];
                if 8*D8 ne ZJ or 4*D8 eq ZJ then
                    continue;
                end if;
                d8ok +:= 1;
                halveTried +:= 1;

                halves := SquareQuarticHalves(g8, v8, f);
                found16 := false;
                for h in halves do
                    G := h[1];
                    ell1q := h[2];
                    for sgn in [1,-1] do
                        D16 := J![G, (sgn*ell1q) mod G];
                        if 2*D16 eq D8 and 16*D16 eq ZJ and 8*D16 ne ZJ then
                            found16 := true;
                            break;
                        end if;
                    end for;
                    if found16 then
                        break;
                    end if;
                end for;
                if not found16 then
                    continue;
                end if;
                halvable +:= 1;

                inv := [];
                try
                    inv := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(fInt))));
                catch e
                    inv := ["torsion_failed"];
                end try;
                contains48 := Type(inv[1]) eq RngIntElt and &or [n mod 48 eq 0 : n in inv];
                print "HALVABLE_AFTER_48_FILTER", "rv", rv, "tv", tv, "bv", bv,
                      "pv", pv, "torsion", inv, "used", used;
                if contains48 then
                    z48Hits +:= 1;
                    print "Z48_TARGET_HIT", "rv", rv, "tv", tv, "bv", bv, "pv", pv, "fInt", fInt;
                end if;
            end for;
        end for;
    end for;
end for;

printf "SEARCH_DONE tested=%o pvals=%o smooth=%o evenSkip=%o factorRich=%o filterPass=%o d8ok=%o halveTried=%o halvable=%o z48Hits=%o\n",
    tested, pvalsFound, smooth, evenSkip, factorRich, filterPass, d8ok, halveTried, halvable, z48Hits;
print "FIRST_KILLS";
for k in Sort([kk : kk in Keys(firstKills)]) do
    printf "  %o : %o\n", k, firstKills[k];
end for;
print "DONE";

