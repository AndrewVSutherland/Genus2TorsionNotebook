
//////////////////////////////////////////////////////////////////////
//  The collision family on the [4,16] cover of M_1(8,4).
//
//  agent_m18_416_branch_discriminant.m found in Res_b(E1core,E0core) a
//  degree-5 factor F5(R,w,a), multiplicity 8, LINEAR in a.  Rational
//  reconstruction from p=101,103,107 gives (cleared by 4):
//
//    4*F5 = 4*(R+1)*(1-w^2)*a
//           + 4*R^5 - 4*R^4 - 8*R^3 + 6*R^2*w^2 + 2*R^2
//           + R*w^4 + 4*R*w^2 - R + 3*w^4 - 6*w^2 - 1 .
//
//  On F5 = 0 the two b-branches of the cover COLLIDE, so the common
//  b-root is rational.  This script:
//   1. verifies symbolically over Q(R,w) that with a = a(R,w) from F5,
//      gcd_b(E1core, E0core) is nontrivial, and extracts the rational
//      double root b(R,w);
//   2. verifies D := f - c4*(x+R)*q^2 = d4*W^2 exactly (W monic quadratic);
//   3. factors d4(R,w) on the family: the remaining condition for a
//      rational half of P_R (an order-16 point) is d4 = SQUARE;
//   4. hunts small rational (R,w) with d4 square and certifies the
//      resulting curves: 2*D16 = +-P_R, torsion, [4,16]/[2,16] check,
//      first-cover tangent condition, and a simplicity certificate.
//
//  Usage: magma -b height:=25 max_hits:=10 agent_m18_416_collision_family.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned height then height := 20;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned max_hits then max_hits := 10;
elif Type(max_hits) eq MonStgElt then max_hits := StringToInteger(max_hits); end if;
if not assigned skip_symbolic then skip_symbolic := false;
elif Type(skip_symbolic) eq MonStgElt then
    skip_symbolic := skip_symbolic in {"true","True","1","yes"};
end if;

//------------------------------------------------------------------
// PART 1: symbolic verification over K2 = Q(R,w)
//------------------------------------------------------------------
K2<R,w> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);
PB<bb> := PolynomialRing(K2);

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*Apol*Bpol;
c4 := R + 2 + 4*t;

aval := (4*R^5 - 4*R^4 - 8*R^3 + 6*R^2*w^2 + 2*R^2
         + R*w^4 + 4*R*w^2 - R + 3*w^4 - 6*w^2 - 1)
        / (4*(R+1)*(w^2-1));

// D as a polynomial in x whose coefficients are polynomials in b (=bb)
PXB<xB> := PolynomialRing(PB);
fB := &+[ PB!Coefficient(f,i) * xB^i : i in [0..5] ];
qB := xB^2 + (PB!aval)*xB + bb;
DB := fB - (PB!c4)*(xB + PB!R)*qB^2;
printf "deg_x D = %o (expect 4)\n", Degree(DB);
d4B := Coefficient(DB,4); d3B := Coefficient(DB,3); d2B := Coefficient(DB,2);
d1B := Coefficient(DB,1); d0B := Coefficient(DB,0);
E1B := 8*d4B^2*d1B - d3B*(4*d4B*d2B - d3B^2);      // in K2[b]
E0B := 64*d4B^3*d0B - (4*d4B*d2B - d3B^2)^2;
printf "deg_b E1 = %o, deg_b E0 = %o\n", Degree(E1B), Degree(E0B);

g := GCD(E1B, E0B);
printf "gcd_b(E1,E0) degree = %o\n", Degree(g);
if Degree(g) lt 1 then
    print "COLLISION FAMILY FAILS: gcd trivial (F5 was an artifact).";
    print "Diagnostic: leading b-coefficients at a=a(R,w):";
    printf "  lc_b E1 = %o\n", LeadingCoefficient(E1B);
    printf "  lc_b E0 = %o\n", LeadingCoefficient(E0B);
    quit;
end if;

rootsB := Roots(g);
printf "rational b-roots of gcd: %o\n", #rootsB;
if #rootsB eq 0 then
    print "gcd nontrivial but no rational root over Q(R,w); gcd =";
    print g;
    quit;
end if;
bval := rootsB[1][1];
printf "b(R,w) = %o\n", bval;

// verify D(bval) = d4*W^2 exactly
q := x^2 + aval*x + bval;
D := f - c4*(x+R)*q^2;
d4v := Coefficient(D,4);
printf "d4 on family = %o\n", d4v;
W := x^2 + (Coefficient(D,3)/(2*d4v))*x
     + (4*d4v*Coefficient(D,2) - Coefficient(D,3)^2)/(8*d4v^2);
issq := D eq d4v*W^2;
printf "D = d4*W^2 exactly: %o\n", issq;
if not issq then
    print "FAMILY VERIFICATION FAILED at the square identity"; quit;
end if;
print "COLLISION FAMILY VERIFIED: rational (a(R,w), b(R,w)) with D = d4*W^2.";

// factor the square condition
d4num := Numerator(d4v); d4den := Denominator(d4v);
print "d4 numerator factorization:";
print Factorization(d4num);
print "d4 denominator factorization:";
print Factorization(d4den);

if skip_symbolic then quit; end if;

//------------------------------------------------------------------
// PART 2: numeric hunt for d4 = square, with full certification
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

function RationalParametersOfHeight(Bd)
    vals := [];
    for den in [1..Bd] do
        for num in [-Bd..Bd] do
            if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

// point-count Frobenius simplicity certificate
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
        if Denominator(dsc) mod pp eq 0 or (Z!Numerator(dsc)) mod pp eq 0 then continue; end if;
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

printf "\nNUMERIC HUNT height=%o\n", height;
params := RationalParametersOfHeight(height);
tested := 0; famOK := 0; d4sq := 0; certified := 0; hits := 0;
torsHist := AssociativeArray();

for Rv in params do
    if Rv in {Q!0, Q!1, Q!-1} then continue; end if;
    for wv in params do
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        tested +:= 1;
        // family data at (Rv,wv)
        den_a := 4*(Rv+1)*(wv^2-1);
        if den_a eq 0 then continue; end if;
        av := (4*Rv^5 - 4*Rv^4 - 8*Rv^3 + 6*Rv^2*wv^2 + 2*Rv^2
               + Rv*wv^4 + 4*Rv*wv^2 - Rv + 3*wv^4 - 6*wv^2 - 1)/den_a;
        bv := Evaluate(Numerator(bval), [Rv,wv]);
        bden := Evaluate(Denominator(bval), [Rv,wv]);
        if bden eq 0 then continue; end if;
        bv := bv/bden;
        // curve data
        tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
        c4v := Rv + 2 + 4*tv;
        if c4v eq 0 then continue; end if;
        fq := xq*(xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4)
                *(c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv));
        if Degree(fq) ne 5 or Discriminant(fq) eq 0 then continue; end if;
        qq := xq^2 + av*xq + bv;
        Dq := fq - c4v*(xq + Rv)*qq^2;
        if Degree(Dq) gt 4 then continue; end if;
        d4q := Coefficient(Dq, 4);
        if d4q eq 0 then continue; end if;
        famOK +:= 1;
        oksq, sq := IsSquareQ(d4q);
        if not oksq then continue; end if;
        d4sq +:= 1;
        // build ell and verify identity
        Wq := xq^2 + (Coefficient(Dq,3)/(2*d4q))*xq
              + (4*d4q*Coefficient(Dq,2) - Coefficient(Dq,3)^2)/(8*d4q^2);
        ellq := sq*Wq;
        if fq - ellq^2 ne c4v*(xq + Rv)*qq^2 then continue; end if;
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
        printf "ORDER16_POINT R=%o w=%o a=%o b=%o d4=%o torsion=%o\n",
            Rv, wv, av, bv, d4q, invs;
        printf "  f = %o\n", fI;
        issimple, pp, chi := SimplicityCertificate(fI);
        printf "  simplicity: %o%o\n",
            issimple select "SIMPLE at p=" cat Sprint(pp) else "no certificate",
            issimple select " chi=" cat Sprint(chi) else "";
        if Has416(invs) then
            hits +:= 1;
            printf "HIT_416 R=%o w=%o torsion=%o simple=%o\n", Rv, wv, invs, issimple;
            if hits ge max_hits then
                printf "STOP max hits\n"; break Rv;
            end if;
        end if;
    end for;
end for;

printf "\nSUMMARY tested=%o famOK=%o d4square=%o certified16=%o hits416=%o\n",
    tested, famOK, d4sq, certified, hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsHist)]) do
    printf "  %o : %o\n", key, torsHist[key];
end for;
print "DONE";
quit;
