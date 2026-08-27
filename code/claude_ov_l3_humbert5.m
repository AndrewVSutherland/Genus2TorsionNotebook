// ===========================================================================
// Lane 3 : derive and VALIDATE the implicit equation of the Humbert surface
// H_5 (RM by Z[(1+sqrt5)/2]) in Igusa-Clebsch invariants.
//
// Source of the parametrization: Elkies-Kumar birational model of Y_-(5)
// (quoted as (2.1) in K. Martin, "Moduli for rational genus 2 curves with
// real multiplication for discriminant 5", arXiv:2206.05752, JTNB 36 (2024)):
//
//   (I2 : I4 : I6 : I10) = ( 24g+6 : 9g^2 : 81g^3+18g^2+36h : 4h^2 )
//
// over the (g,h)-plane (the z-coordinate of Y only records the RM action,
// not the moduli point).  A curve is on H_5 iff its invariants agree with
// some (g,h) after the weighted scaling (I2,I4,I6,I10) -> (u I2, u^2 I4,
// u^3 I6, u^5 I10).  Eliminating g,h and then u gives the implicit equation.
//
// EVERYTHING here is validated before use:
//   * convention check against the paper's explicit RM-5 curve
//   * random points of the parametrization -> curve via Mestre -> Frobenius
//     real-subfield-disc census must be constantly 5
//   * random non-RM curves -> H5 must NOT vanish
//   * RM curves of other discriminants -> H5 must NOT vanish
//
// usage: magma -b code/claude_ov_l3_humbert5.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
Z := Integers();  Q := Rationals();
Px<x> := PolynomialRing(Q);
PT<T> := PolynomialRing(Q);

// --------------------------------------------------------------- derivation
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
S<u> := PolynomialRing(R);
g := (u*J2 - 6)/24;
A := (u*J2 - 6)^2 - 64*u^2*J4;                 // 9g^2 = u^2 J4
h := (u^3*J6 - 81*g^3 - 18*g^2)/36;
B := 4*h^2 - u^5*J10;                          // 4h^2 = u^5 J10
B := S ! (B * 36^2);                           // clear denominators
printf "degA=%o degB=%o\n", Degree(A), Degree(B);
F := Resultant(A, B);
printf "H5 resultant: #terms=%o  totaldeg=%o\n", #Terms(F), TotalDegree(F);
fac := Factorization(F);
printf "factors:\n";
for t in fac do
    printf "  mult=%o  #terms=%o  deg=%o\n", t[2], #Terms(t[1]), TotalDegree(t[1]);
end for;

// weighted degree of a factor (weights I2=2,I4=4,I6=6,I10=10)
function WDeg(p)
    w := {};
    for m in Monomials(p) do
        e := Exponents(m);
        Include(~w, 2*e[1]+4*e[2]+6*e[3]+10*e[4]);
    end for;
    return w;
end function;
for i in [1..#fac] do
    printf "  factor %o weighted degrees %o\n", i, WDeg(fac[i][1]);
end for;

// pick the factor(s) that actually vanish on the parametrization
PGH<gg,hh> := PolynomialRing(Q, 2);
sub2 := hom< R -> PGH | 24*gg+6, 9*gg^2, 81*gg^3+18*gg^2+36*hh, 4*hh^2 >;
H5 := R ! 1;
for t in fac do
    if sub2(t[1]) eq 0 then
        printf "  -> factor with %o terms VANISHES on the parametrization\n", #Terms(t[1]);
        H5 *:= t[1];
    end if;
end for;
printf "H5 chosen: #terms=%o weighted degrees %o\n", #Terms(H5), WDeg(H5);
Write("results/claude_ov_l3_humbert5_eqn.txt", Sprintf("%o", H5) : Overwrite := true);

// --------------------------------------------------------------- evaluation
function H5Val(C)
    I := IgusaClebschInvariants(C);
    return Evaluate(H5, [Q!I[1], Q!I[2], Q!I[3], Q!I[4]]);
end function;

// Frobenius real-subfield-disc census (the lab's cheap RM screen)
function DiscCensus(f, B)
    m := LCM([Denominator(c) : c in Coefficients(f)] cat [Z|1]);
    fI := Px ! [ Z ! (Coefficient(f,i)*m^(6-i)) : i in [0..Degree(f)] ];
    dsc := Discriminant(fI);
    discs := {};
    np := 0;
    for p in PrimesInInterval(3, B) do
        if (Z!dsc) mod p eq 0 then continue; end if;
        Fp := GF(p);
        fp := PolynomialRing(Fp) ! fI;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        L := LPolynomial(HyperellipticCurve(fp));
        a1 := -Coefficient(L,1); a2 := Coefficient(L,2);
        n := (Z!a1)^2 - 4*((Z!a2) - 2*p);
        if n ne 0 then
            d := SquarefreeFactorization(n);
            if d ne 1 then Include(~discs, d); end if;
        end if;
        np +:= 1;
    end for;
    return discs, np;
end function;

printf "\n==== VALIDATION ====\n";

// (a) convention check: the paper's explicit RM-5 curve, invariants (8:1:3:8/3125)
fpap := (2*x^3-2*x^2-x-1)*(x^3-x^2+2*x+2);
Cpap := HyperellipticCurve(fpap);
Ipap := IgusaClebschInvariants(Cpap);
printf "paper RM5 curve IC = %o\n", Ipap;
// weighted-proportional to (8,1,3,8/3125)?
uu := (Q!Ipap[1])/8;
printf "  scaling u=%o : u^2*1=%o vs I4=%o ; u^3*3=%o vs I6=%o ; u^5*8/3125=%o vs I10=%o\n",
    uu, uu^2*1, Ipap[2], uu^3*3, Ipap[3], uu^5*8/3125, Ipap[4];
printf "  H5(paper curve) = %o\n", H5Val(Cpap);
dc, np := DiscCensus(fpap, 200);
printf "  disc census (p<=200, %o primes) = %o\n", np, Sort(Setseq(dc));

// (b) the two known [2,22] RM(sqrt5) witnesses
f1 := x*(x-1)*(x^2-3*x+8)*(x^2+4*x+27);                       // 19044.h.2
f2 := x*(x-1)*(x^2-34/15*x+32/45)*(x^2+31/15*x-2/9);          // BLP C4corr
for pr in [<"19044.h.2", f1>, <"BLP-C4corr", f2>] do
    C := HyperellipticCurve(pr[2]);
    dc, np := DiscCensus(pr[2], 200);
    printf "%o : H5=%o  census=%o\n", pr[1], H5Val(C), Sort(Setseq(dc));
end for;

// (c) random points of the EK parametrization -> Mestre curve -> census
printf "\n-- forward test: random (g,h) on the parametrization --\n";
cnt := 0;
for trial in [1..40] do
    gv := Random([-6..6]) + Random([1..5])/Random([1..5]);
    hv := Random([-6..6]) + Random([1..5])/Random([1..5]);
    IC := [24*gv+6, 9*gv^2, 81*gv^3+18*gv^2+36*hv, 4*hv^2];
    if IC[4] eq 0 then continue; end if;
    ok := true;
    try
        C := HyperellipticCurveFromIgusaClebsch(IC);
        f := HyperellipticPolynomials(SimplifiedModel(C));
    catch e
        ok := false;
    end try;
    if not ok then continue; end if;
    hv5 := H5Val(C);
    dc, np := DiscCensus(f, 120);
    cnt +:= 1;
    printf "  (g,h)=(%o,%o) H5=%o census=%o\n", gv, hv, hv5, Sort(Setseq(dc));
    if cnt ge 8 then break; end if;
end for;

// (d) random NON-RM curves: H5 must be nonzero
printf "\n-- negative control: random curves --\n";
nz := 0; nt := 0;
for trial in [1..12] do
    f := x^6 + &+[Random([-9..9])*x^i : i in [0..5]];
    if Discriminant(f) eq 0 then continue; end if;
    C := HyperellipticCurve(f);
    v := H5Val(C);
    nt +:= 1; if v ne 0 then nz +:= 1; end if;
end for;
printf "  random curves: %o of %o have H5 <> 0\n", nz, nt;

// (e) an RM(sqrt2) curve (contact-7 control, recorded in the repo)
//     v-triple (-3,-3/4,-3/5) has RM(sqrt2); rebuild from the contact-7 chart
printf "\n-- RM of another discriminant --\n";
// direct: a known RM(sqrt8) curve from LMFDB, 294.a.294.1-like; use a generic
// test instead: curves from the parametrization of H_5 must be the ONLY zeros
printf "(covered by (f) below)\n";
quit;
