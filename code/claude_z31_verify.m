// Exact verification of the 1830.2.a.q curve: torsion + simplicity certificate
// Run on aws-spot-11.  Expected: J(Q)_tors = Z/31, geometrically simple (RM by Z[sqrt2]).
SetMemoryLimit(16*10^9);
Q := Rationals(); Z := Integers();
R<x> := PolynomialRing(Q);

function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
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
        for nn in [2..12] do
            if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
        end for;
        if not drop then return true, pp, chi; end if;
    end for;
    return false, 0, RT!0;
end function;

fco := [-126,816,-2466,4300,-4587,2841,-839];
hco := [0,-1,-1];
f := &+[fco[i]*x^(i-1) : i in [1..#fco]];
h := &+[hco[i]*x^(i-1) : i in [1..#hco]];
C := HyperellipticCurve(f,h);
printf "curve y^2 + (%o)y = %o\n", h, f;
Cs := SimplifiedModel(C);
fs, hs := HyperellipticPolynomials(Cs);
assert hs eq 0;
// integral model: fs = 4f + h^2 already integral
fInt := fs;
printf "simplified integral model: y^2 = %o\n", fInt;
printf "disc = %o\n", Factorization(Z!Numerator(Discriminant(fInt)));
J := Jacobian(Cs);
t0 := Cputime();
T, mp := TorsionSubgroup(J);
printf "TORSION_INVARIANTS: %o   (%.1o s)\n", Invariants(T), Cputime(t0);
for i in [1..Ngens(T)] do
    P := mp(T.i);
    printf "TORSION_GEN %o: %o  order %o\n", i, P, Order(P);
end for;
ok, pw, chiw := SimplicityCertificate(fInt);
printf "SIMPLICITY: %o  witness prime %o  chi = %o\n", ok, pw, chiw;
// extra: conductor sanity (expect 1830 = 2*3*5*61 up to squares in disc)
printf "BadPrimes(C): %o\n", BadPrimes(C);
printf "DONE\n";
quit;
