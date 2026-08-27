
//////////////////////////////////////////////////////////////////////
//  Main-sheet analysis for [4,16], part 1: the a = 0 slice of the
//  degree-52 hypersurface F52(R,w,a) = 0.
//
//  Background: Res_b(E1core, E0core) = (w^2-1)^3 (R+1)^8 F5^8 F52 with
//  F5 the degenerate-branch (d4=0) locus.  The main sheet of the cover
//  Sigma' projects birationally onto F52 = 0 (b = root of the linear
//  gcd, hence RATIONAL for free); the remaining lift condition is
//  d4 = square.  The <0,0> 7-adic tower lives on the main sheet with
//  a == 0 mod 7, and the origin (R,w,a,b)=(0,0,0,0) is on the main
//  sheet only (d4(origin)=1).
//
//  This script computes the a=0 slice exactly over Q:
//    Res_b(E1|_{a=0}, E0|_{a=0}) in Q[R,w],
//  strips the known boundary/F5-slice factors, factors the remainder
//  (= F52(R,w,0)), reports genus/rationality of components -- in
//  particular the component(s) through the origin (the tower shadow) --
//  and hunts rational points with the full order-16 certification.
//
//  Also reports the degree profile of F52 and its mod-7 factorization
//  (tower shadow) from a mod-p resultant.
//
//  Usage: magma -b height:=30 agent_m18_416_f52_a0_slice.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned height then height := 25;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned do_modp then do_modp := true;
elif Type(do_modp) eq MonStgElt then do_modp := do_modp in {"true","True","1","yes"}; end if;

// ---- build E1, E0 in Q[R,w,a,b] ----
R4<R,w,a,b> := PolynomialRing(Q, 4, "grevlex");
K4 := FieldOfFractions(R4);
PX<x> := PolynomialRing(K4);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*Apol*Bpol;
c4 := R + 2 + 4*t;
q := x^2 + a*x + b;
D := f - c4*(x+R)*q^2;
d4 := Coefficient(D,4); d3 := Coefficient(D,3); d2 := Coefficient(D,2);
d1 := Coefficient(D,1); d0 := Coefficient(D,0);
E1 := 8*d4^2*d1 - d3*(4*d4*d2 - d3^2);
E0 := 64*d4^3*d0 - (4*d4*d2 - d3^2)^2;
function ClearPrim(g)
    gn := Numerator(g);
    den := LCM([Denominator(co) : co in Coefficients(gn)]);
    gn := den*gn;
    gg := GCD([Z!co : co in Coefficients(gn)]);
    return R4!(gn/gg);
end function;
E1n := ClearPrim(E1); E0n := ClearPrim(E0);
while IsDivisibleBy(E1n, R-1) do E1n := ExactQuotient(E1n, R-1); end while;
while IsDivisibleBy(E0n, R-1) do E0n := ExactQuotient(E0n, R-1); end while;

// F5 = numerator(d4)/(R-1) for slice stripping
d4n := ClearPrim(d4);
F5 := ExactQuotient(d4n, R-1);

// ---- degree profile + mod-7 shadow of F52 (mod-p resultants) ----
if do_modp then
    for pp in [101, 7] do
        Fp := GF(pp);
        R3p<Rp,wp,ap> := PolynomialRing(Fp, 3, "grevlex");
        Pb := PolynomialRing(R3p);
        hp := hom<R4 -> Pb | [Pb!Rp, Pb!wp, Pb!ap, Pb.1]>;
        E1p := hp(E1n); E0p := hp(E0n);
        Res := Resultant(E1p, E0p);
        if Res eq 0 then printf "p=%o: resultant vanishes identically\n", pp; continue; end if;
        facR := Factorization(Res);
        printf "p=%o Res_b factor profile:\n", pp;
        for fr in facR do
            printf "  <deg %o, degR %o, degw %o, dega %o, mult %o>",
                TotalDegree(fr[1]), Degree(fr[1],1), Degree(fr[1],2),
                Degree(fr[1],3), fr[2];
            if #Terms(fr[1]) le 8 then printf "  %o", fr[1]; end if;
            printf "\n";
        end for;
    end for;
end if;

// ---- exact a=0 slice ----
print "computing exact a=0 slice resultant over Q ...";
E1s := Evaluate(E1n, [R, w, R4!0, b]);
E0s := Evaluate(E0n, [R, w, R4!0, b]);
tim := Cputime();
ResS := Resultant(E1s, E0s, b);
printf "Res_b on a=0 slice: total degree %o, terms %o (%o s)\n",
    TotalDegree(ResS), #Terms(ResS), Cputime(tim);

// strip known factors: (w^2-1) powers, (R+1) powers, F5(a=0) powers
F5s := Evaluate(F5, [R, w, R4!0, b]);
for gfac in [R4!(w-1), R4!(w+1), R4!(R+1), R4!R, R4!w, F5s] do
    while IsDivisibleBy(ResS, gfac) do ResS := ExactQuotient(ResS, gfac); end while;
end for;
printf "after stripping boundary/F5 factors: degree %o, terms %o\n",
    TotalDegree(ResS), #Terms(ResS);

print "factorization of the a=0 slice core (F52(R,w,0) up to strips):";
facS := Factorization(ResS);
for ff in facS do
    g := ff[1];
    printf "  <deg %o, degR %o, degw %o, terms %o, mult %o> origin=%o",
        TotalDegree(g), Degree(g,1), Degree(g,2), #Terms(g), ff[2],
        Evaluate(g, [0,0,0,0]) eq 0;
    if #Terms(g) le 10 then printf "  %o", g; end if;
    printf "\n";
end for;

// genus of components (projective closure), especially through origin
A2<Ru, wu> := AffineSpace(Q, 2);
for ff in facS do
    g := ff[1];
    dg := TotalDegree(g);
    if dg le 1 then printf "component deg 1: rational line %o\n", g; continue; end if;
    if dg gt 32 then printf "component deg %o: genus skipped (too big)\n", dg; continue; end if;
    cf := Evaluate(g, [Ru, wu, 0, 0]);
    Cv := Curve(A2, cf);
    try
        gg := Genus(ProjectiveClosure(Cv));
        printf "component deg %o (origin=%o): genus %o\n",
            dg, Evaluate(g,[0,0,0,0]) eq 0, gg;
    catch err
        printf "component deg %o: genus failed (%o)\n", dg, err`Object;
    end try;
end for;
print "DONE_SYMBOLIC";

// ---- hunt rational points on the slice core, certify order 16 ----
P<xq> := PolynomialRing(Q);
Pw<wp2> := PolynomialRing(Q);

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
    for den in [1..Bd] do for num in [-Bd..Bd] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

printf "\nHUNT on a=0 slice, height=%o\n", height;
params := RationalParametersOfHeight(height);
swept := 0; pts := 0; certified := 0; hits := 0;
torsHist := AssociativeArray();

for Rv in params do
    if Rv in {Q!0, Q!1, Q!-1} then continue; end if;
    swept +:= 1;
    // specialize slice core at R=Rv: polynomial in w
    Sw := Pw!0;
    for tm in Terms(ResS) do
        dR := Degree(tm, R); dW := Degree(tm, w);
        Sw +:= LeadingCoefficient(tm)*Rv^dR*wp2^dW;
    end for;
    if Sw eq 0 then continue; end if;
    for rt in Roots(Sw) do
        wv := rt[1];
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        pts +:= 1;
        // main-sheet point at (Rv, wv, a=0): find rational b via gcd
        tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
        c4v := Rv + 2 + 4*tv;
        if c4v eq 0 then continue; end if;
        fq := xq*(xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4)
                *(c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv));
        if Degree(fq) ne 5 or Discriminant(fq) eq 0 then continue; end if;
        // b-polys at (Rv,wv,a=0)
        Pb2<bq> := PolynomialRing(Q);
        E1v := Pb2!0; E0v := Pb2!0;
        for tm in Terms(E1s) do
            E1v +:= LeadingCoefficient(tm)*Rv^Degree(tm,R)*wv^Degree(tm,w)*bq^Degree(tm,b);
        end for;
        for tm in Terms(E0s) do
            E0v +:= LeadingCoefficient(tm)*Rv^Degree(tm,R)*wv^Degree(tm,w)*bq^Degree(tm,b);
        end for;
        gv := GCD(E1v, E0v);
        if Degree(gv) lt 1 then continue; end if;
        for rtb in Roots(gv) do
            bv := rtb[1];
            qq := xq^2 + bv;   // a = 0
            Dq := fq - c4v*(xq + Rv)*qq^2;
            if Degree(Dq) ne 4 then continue; end if;
            d4q := Coefficient(Dq, 4);
            if d4q eq 0 then continue; end if;
            // square quartic check + lift condition
            Wq := xq^2 + (Coefficient(Dq,3)/(2*d4q))*xq
                  + (4*d4q*Coefficient(Dq,2) - Coefficient(Dq,3)^2)/(8*d4q^2);
            if Dq ne d4q*Wq^2 then continue; end if;
            oksq, sq := IsSquareQ(d4q);
            printf "MAIN_SHEET_POINT R=%o w=%o b=%o d4=%o d4_square=%o\n",
                Rv, wv, bv, d4q, oksq;
            if not oksq then continue; end if;
            ellq := sq*Wq;
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
end for;
printf "\nSUMMARY swept_R=%o slice_points=%o certified16=%o hits416=%o\n",
    swept, pts, certified, hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsHist)]) do
    printf "  %o : %o\n", key, torsHist[key];
end for;
print "DONE";
quit;
