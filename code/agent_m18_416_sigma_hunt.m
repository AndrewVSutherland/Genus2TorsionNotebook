
//////////////////////////////////////////////////////////////////////
//  Rational-point hunt on the degenerate-branch curve Sigma(R,w)=0 of
//  the [4,16] cover (see agent_m18_416_degenerate_branch.m), skipping
//  the (infeasible) genus computation of the degree-40 core.
//
//  Sigma is even in w, so substitute W=w^2: for each rational R, find
//  rational roots W of a half-degree polynomial, keep W = square, then
//  run the full order-16 certification chain.
//
//  Also reports mod-p point counts of the core curve (plausibility).
//
//  Usage: magma -b height:=60 agent_m18_416_sigma_hunt.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned height then height := 40;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned modp_counts then modp_counts := true;
elif Type(modp_counts) eq MonStgElt then
    modp_counts := modp_counts in {"true","True","1","yes"};
end if;

// ---- symbolic branch data over Q(R,w) ----
K2<R,w> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*Apol*Bpol;
c4 := R + 2 + 4*t;
f4 := Coefficient(f, 4);
aval := (f4/c4 - R)/2;
f3 := Coefficient(f, 3);
bval := (f3/c4 - aval^2 - 2*aval*R)/2;
q := x^2 + aval*x + bval;
D := f - c4*(x+R)*q^2;
d2 := Coefficient(D,2); d1 := Coefficient(D,1); d0 := Coefficient(D,0);
Sigma := d1^2 - 4*d2*d0;
Snum := Numerator(Sigma);
// strip R-1 boundary factors
RW := Parent(Snum);
Rv_ := RW.1; wv_ := RW.2;
while IsDivisibleBy(Snum, Rv_ - 1) do Snum := ExactQuotient(Snum, Rv_ - 1); end while;
printf "Sigma core: degree %o, degR %o, degw %o, terms %o\n",
    TotalDegree(Snum), Degree(Snum,1), Degree(Snum,2), #Terms(Snum);

// check evenness in w
isEvenW := &and[ IsEven(Degree(tm, wv_)) : tm in Terms(Snum) ];
printf "Sigma even in w: %o\n", isEvenW;
error if not isEvenW, "expected even structure in w";

// mod-p plausibility counts of the affine core curve
if modp_counts then
    for pp in [11,13,17,19,23] do
        Fp := GF(pp);
        cnt := 0; cntd2sq := 0;
        for r0 in Fp do for w0 in Fp do
            if Evaluate(Snum, [r0,w0]) ne 0 then continue; end if;
            if r0 in {Fp!0,Fp!1,Fp!-1} or w0 in {Fp!0,Fp!1,Fp!-1} then continue; end if;
            cnt +:= 1;
            dd := Denominator(d2);
            if Evaluate(dd, [r0,w0]) eq 0 then continue; end if;
            nn := Evaluate(Numerator(d2), [r0,w0]);
            if nn eq 0 then continue; end if;
            if IsSquare(nn/Evaluate(dd,[r0,w0])) then cntd2sq +:= 1; end if;
        end for; end for;
        printf "p=%o: affine core points %o, with d2 square %o\n", pp, cnt, cntd2sq;
    end for;
end if;

// ---- hunt: sweep R, solve for W=w^2 ----
P<xq> := PolynomialRing(Q);
Pw<Wv> := PolynomialRing(Q);

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false, 0; end if;
    okn, sn := IsSquare(Numerator(qv));
    okd, sd := IsSquare(Denominator(qv));
    if okn and okd then return true, sn/sd; end if;
    return false, 0;
end function;

function IntegralModelPolynomial(fq)
    L := 1;
    for i in [0..Degree(fq)] do L := LCM(L, Denominator(Coefficient(fq, i))); end for;
    return P!(L^2*fq), L;
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
    for den in [1..Bd] do
        for num in [-Bd..Bd] do
            if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

printf "\nHUNT height=%o (sweep R, roots in W=w^2, W square)\n", height;
params := RationalParametersOfHeight(height);
swept := 0; Wroots := 0; wRat := 0; certified := 0; hits := 0;
torsHist := AssociativeArray();

for Rv0 in params do
    if Rv0 in {Q!0, Q!1, Q!-1} then continue; end if;
    swept +:= 1;
    // specialize Snum at R=Rv0 as a polynomial in W = w^2
    Sw := Pw!0;
    for tm in Terms(Snum) do
        dR := Degree(tm, Rv_); dW := Degree(tm, wv_) div 2;
        Sw +:= LeadingCoefficient(tm)*Rv0^dR*Wv^dW;
    end for;
    if Sw eq 0 then continue; end if;
    for rt in Roots(Sw) do
        Wq := rt[1];
        Wroots +:= 1;
        okW, wv0 := IsSquareQ(Wq);
        if not okW then continue; end if;
        if wv0 in {Q!0, Q!1} then continue; end if;
        wRat +:= 1;
        // certification at (Rv0, wv0)
        Rv := Rv0; wv := wv0;
        tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
        c4v := Rv + 2 + 4*tv;
        if c4v eq 0 then continue; end if;
        fq := xq*(xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4)
                *(c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv));
        if Degree(fq) ne 5 or Discriminant(fq) eq 0 then continue; end if;
        f4v := Coefficient(fq, 4);
        av := (f4v/c4v - Rv)/2;
        f3v := Coefficient(fq, 3);
        bv := (f3v/c4v - av^2 - 2*av*Rv)/2;
        qq := xq^2 + av*xq + bv;
        Dq := fq - c4v*(xq + Rv)*qq^2;
        if Degree(Dq) gt 2 then continue; end if;
        d2q := Coefficient(Dq, 2);
        if d2q eq 0 then continue; end if;
        oksq, sq := IsSquareQ(d2q);
        printf "SIGMA_POINT R=%o w=%o d2=%o d2_square=%o\n", Rv, wv, d2q, oksq;
        if not oksq then continue; end if;
        ellq := sq*(xq + Coefficient(Dq,1)/(2*d2q));
        if fq - ellq^2 ne c4v*(xq + Rv)*qq^2 then continue; end if;
        fI, L := IntegralModelPolynomial(fq);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        Qfac := Rv^2 - (Q!1/2)*Rv*wv^2 + (Q!1/2)*Rv - wv^2;
        YR := -2*Rv*(Rv-1)^2*Qfac/(wv^2-1);
        if Evaluate(fI, -Rv) ne (L*YR)^2 then continue; end if;
        PR := J![xq + Rv, P!(L*YR)];
        vred := (L*ellq) mod qq;
        good := false; D16 := J!0;
        for sgn in [1,-1] do
            try D16c := J![qq, sgn*vred];
            catch err continue; end try;
            if 2*D16c eq PR or 2*D16c eq -PR then good := true; D16 := D16c; break; end if;
        end for;
        if not good then printf "  halving check failed\n"; continue; end if;
        if 16*D16 ne J!0 or 8*D16 eq J!0 then printf "  order not 16\n"; continue; end if;
        certified +:= 1;
        invs := Invariants(TorsionSubgroup(J));
        key := Sprint(invs);
        if IsDefined(torsHist,key) then torsHist[key] +:= 1; else torsHist[key] := 1; end if;
        printf "ORDER16_POINT R=%o w=%o torsion=%o\n", Rv, wv, invs;
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

printf "\nSUMMARY swept_R=%o Wroots=%o w_rational=%o certified16=%o hits416=%o\n",
    swept, Wroots, wRat, certified, hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsHist)]) do
    printf "  %o : %o\n", key, torsHist[key];
end for;
print "DONE";
quit;
