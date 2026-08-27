
//////////////////////////////////////////////////////////////////////
//  The degenerate-ell branch of the [4,16] cover of M_1(8,4), reduced
//  to a plane curve.
//
//  For the halving f - ell^2 = c4*(x+R)*q^2 of P_R = (-R,Y_R), set
//  D := f - c4*(x+R)*q^2 (deg <= 4).  The branch where ell drops to
//  degree 1 is:
//      d4 = 0  -> solves a = a(R,w)        (d4 linear in a),
//      d3 = 0  -> solves b = b(R,w)        (d3 linear in b),
//  leaving D = d2*x^2 + d1*x + d0 and the halving condition
//      D = ell^2, deg(ell) = 1  <=>  Sigma := d1^2 - 4*d2*d0 = 0
//                                    and d2 a rational SQUARE.
//  Sigma(R,w) = 0 is an explicit PLANE CURVE: every rational point on
//  it with d2 square gives a rational half of P_R, i.e. a rational
//  order-16 point on the corresponding genus-2 Jacobian in M_1(8,4)
//  (torsion at least [2,16]; with the first tangent cover, [4,16]).
//
//  This script derives Sigma, factors it, reports the genus/rationality
//  of its components, and hunts+certifies rational points.
//
//  Usage: magma -b height:=40 agent_m18_416_degenerate_branch.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned height then height := 30;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned max_hits then max_hits := 20;
elif Type(max_hits) eq MonStgElt then max_hits := StringToInteger(max_hits); end if;

//------------------------------------------------------------------
// symbolic derivation over Q(R,w)
//------------------------------------------------------------------
K2<R,w> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*Apol*Bpol;
c4 := R + 2 + 4*t;

// generic q with symbolic a,b: build D coefficients directly
// d4 = f4 - c4*(2a + R), so a(R,w) from d4=0:
f4 := Coefficient(f, 4);
aval := (f4/c4 - R)/2;
// with a = aval: d3 = f3 - c4*( (a^2 + 2b) + 2aR ), so b(R,w) from d3=0:
f3 := Coefficient(f, 3);
bval := ((f3/c4 - aval^2 - 2*aval*R))/2;
printf "a(R,w) = %o\n\n", aval;
printf "b(R,w) = %o\n\n", bval;

q := x^2 + aval*x + bval;
D := f - c4*(x+R)*q^2;
printf "deg D on branch = %o (expect 2)\n", Degree(D);
d2 := Coefficient(D,2); d1 := Coefficient(D,1); d0 := Coefficient(D,0);

Sigma := d1^2 - 4*d2*d0;
Snum := Numerator(Sigma);
Sden := Denominator(Sigma);
printf "Sigma numerator: total degree %o, terms %o\n",
    TotalDegree(Snum), #Terms(Snum);
print "Sigma numerator factorization:";
facS := Factorization(Snum);
for ff in facS do
    printf "  <deg %o, degR %o, degw %o, terms %o, mult %o>",
        TotalDegree(ff[1]), Degree(ff[1],1), Degree(ff[1],2), #Terms(ff[1]), ff[2];
    if #Terms(ff[1]) le 12 then printf "  %o", ff[1]; end if;
    printf "\n";
end for;

print "d2 numerator factorization:";
d2num := Numerator(d2); d2den := Denominator(d2);
print Factorization(d2num);
print "d2 denominator factorization:";
print Factorization(d2den);

// genus / rationality of the interesting Sigma components
print "component analysis:";
A2<Ra, wa> := AffineSpace(Q, 2);
for ff in facS do
    g := ff[1];
    dg := TotalDegree(g);
    if dg le 1 then
        printf "  deg-1 component %o : rational line\n", g;
        continue;
    end if;
    // map into the affine plane
    cf := Evaluate(g, [A2.1, A2.2]);
    Cv := Curve(A2, cf);
    try
        gg := Genus(ProjectiveClosure(Cv));
        printf "  deg-%o component: genus %o\n", dg, gg;
    catch err
        printf "  deg-%o component: genus computation failed (%o)\n", dg, err`Object;
    end try;
end for;

//------------------------------------------------------------------
// rational point hunt on Sigma = 0 with d2 = square, plus certification
//------------------------------------------------------------------
P<xq> := PolynomialRing(Q);

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

// enumerate rational points on Sigma=0: sweep R, solve numerator for w
printf "\nHUNT height=%o: sweep R, solve Sigma(R,.)=0 for rational w\n", height;
RW3 := PolynomialRing(Q, 2);
SnumP := Evaluate(Snum, [RW3.1, RW3.2]);
Pw<wp> := PolynomialRing(Q);

function RationalParametersOfHeight(Bd)
    vals := [];
    for den in [1..Bd] do
        for num in [-Bd..Bd] do
            if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

params := RationalParametersOfHeight(height);
swept := 0; ptsFound := 0; d2sq := 0; certified := 0; hits := 0;
torsHist := AssociativeArray();
seen := {};

for Rv in params do
    if Rv in {Q!0, Q!1, Q!-1} then continue; end if;
    swept +:= 1;
    // specialize Sigma numerator at R = Rv: polynomial in w
    Sw := Pw!0;
    for tm in Terms(SnumP) do
        dR := Degree(tm, RW3.1); dW := Degree(tm, RW3.2);
        Sw +:= LeadingCoefficient(tm)*Rv^dR*wp^dW;
    end for;
    if Sw eq 0 then continue; end if;
    rts := Roots(Sw);
    for rt in rts do
        wv := rt[1];
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        if <Rv,wv> in seen then continue; end if;
        Include(~seen, <Rv,wv>);
        ptsFound +:= 1;
        // family data
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
        // Sigma=0 should make Dq a perfect square up to d2:
        oksq, sq := IsSquareQ(d2q);
        if not oksq then continue; end if;
        d2sq +:= 1;
        ellq := sq*(xq + Coefficient(Dq,1)/(2*d2q));
        if fq - ellq^2 ne c4v*(xq + Rv)*qq^2 then
            printf "IDENTITY FAIL R=%o w=%o\n", Rv, wv; continue;
        end if;
        // Jacobian certification
        fI, L := IntegralModelPolynomial(fq);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        Qfac := Rv^2 - (Q!1/2)*Rv*wv^2 + (Q!1/2)*Rv - wv^2;
        YR := -2*Rv*(Rv-1)^2*Qfac/(wv^2-1);
        if Evaluate(fI, -Rv) ne (L*YR)^2 then continue; end if;
        PR := J![xq + Rv, P!(L*YR)];
        ellI := L*ellq;
        vred := ellI mod qq;
        good := false; D16 := J!0;
        for sgn in [1,-1] do
            try
                D16c := J![qq, sgn*vred];
            catch err continue; end try;
            if 2*D16c eq PR or 2*D16c eq -PR then
                good := true; D16 := D16c; break;
            end if;
        end for;
        if not good then continue; end if;
        if 16*D16 ne J!0 or 8*D16 eq J!0 then continue; end if;
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
        if hits ge max_hits then break Rv; end if;
    end for;
end for;

printf "\nSUMMARY swept_R=%o Sigma_points=%o d2square=%o certified16=%o hits416=%o\n",
    swept, ptsFound, d2sq, certified, hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsHist)]) do
    printf "  %o : %o\n", key, torsHist[key];
end for;
print "DONE";
quit;
