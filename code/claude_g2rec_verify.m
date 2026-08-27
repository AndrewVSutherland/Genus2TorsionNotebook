// claude_g2rec_verify.m -- full verification of a reconstructed modular genus-2 curve:
//   (1) L-function match: Euler factors of Jac(C) vs exact newform a_p data
//       (APDATA lines in the mdeep log) at every good p <= 199, with automatic
//       quadratic-twist detection if the raw curve mismatches;
//   (2) exact TorsionSubgroup of the (correctly twisted) Jacobian;
//   (3) geometric-simplicity certificate (D4 / root-power criterion, claude_z31_verify.m);
//   (4) repaired RM screen (squarefree core of a1^2-4(a2-2p), a1<>0, non-square)
//       -- expect {5} for RM by Q(sqrt5).
// Usage: magma -b Lab:=... N:=... FC:=c0,c1,... HCC:=c0,... APF:=mdeep_XXXX.log claude_g2rec_verify.m
SetColumns(0);
SetMemoryLimit(60*10^9);
Q := Rationals(); Z := Integers();
R<x> := PolynomialRing(Q);
N := StringToInteger(N);
fco := [StringToInteger(s) : s in Split(FC, ",")];
f := &+[fco[i]*x^(i-1) : i in [1..#fco]];
if assigned HCC and HCC ne "0" then
    hco := [StringToInteger(s) : s in Split(HCC, ",")];
    h := &+[hco[i]*x^(i-1) : i in [1..#hco]];
else
    h := R!0;
end if;
printf "VERIFY %o : y^2 + (%o)y = %o\n", Lab, h, f;
C := HyperellipticCurve(f, h);
Cs := SimplifiedModel(C);
fs, hs := HyperellipticPolynomials(Cs);
assert hs eq 0;
// clear denominators to an integral model y^2 = fInt
den := LCM([Denominator(c) : c in Coefficients(fs)]);
fInt := fs * den^2;
CInt := HyperellipticCurve(fInt);
printf "integral model: y^2 = %o\n", fInt;

// ---- parse APDATA ----
ap := AssociativeArray();
for L in Split(Read(APF), "\n") do
    if #L ge 7 and L[1..6] eq "APDATA" then
        s := Split(L, " ");
        ap[StringToInteger(s[2])] := <StringToInteger(s[3]), StringToInteger(s[4])>;
    end if;
end for;
plist := Sort(Setseq(Keys(ap)));
printf "APDATA primes: %o (n=%o)\n", plist, #plist;

DC := Z!Discriminant(fInt);
// Euler factor of the curve at good p, as polynomial in x over Z
function LpolyC(fI, p)
    Jp := Jacobian(ChangeRing(HyperellipticCurve(fI), GF(p)));
    return R!EulerFactor(Jp);
end function;
function TargetLpoly(p, tp, np)
    return 1 - tp*x + (np + 2*p)*x^2 - p*tp*x^3 + p^2*x^4;
end function;

// ---- twist detection: find squarefree d with a_p(C^d) = chi_d(p) a_p(f) ----
cands := [1];
badps := [pr[1] : pr in Factorization(2*N*DC)];
// all squarefree products of bad primes, both signs, capped
prods := [1];
for pq in badps do
    prods := prods cat [Z!(t*pq) : t in prods];
    if #prods gt 512 then break; end if;
end for;
cands := Sort(Setseq(Seqset(prods cat [-t : t in prods])));
bestd := 0; bestgood := -1;
for d in cands do
    good := 0; ok := true;
    for p in plist do
        if DC mod p eq 0 or d mod p eq 0 then continue; end if;
        tp, np := Explode(ap[p]);
        chi := KroneckerSymbol(d, p);
        Ep := LpolyC(fInt, p);
        if Ep eq TargetLpoly(p, chi*tp, np) then
            good +:= 1;
        else
            ok := false; break;
        end if;
    end for;
    if ok and good gt bestgood then bestgood := good; bestd := d; end if;
    if ok and d eq 1 then break; end if; // identity twist works: done
end for;
if bestd eq 0 then
    printf "LMATCH_FAIL %o : no quadratic twist matches all good primes\n", Lab;
    quit;
end if;
printf "TWIST d=%o matches at %o good primes (of %o APDATA primes)\n", bestd, bestgood, #plist;
Cfin := bestd eq 1 select CInt else QuadraticTwist(CInt, bestd);
Cfin := ReducedMinimalWeierstrassModel(Cfin);
ffin, hfin := HyperellipticPolynomials(SimplifiedModel(Cfin));
fver := 4*ffin + hfin^2;
denv := LCM([Denominator(c) : c in Coefficients(fver)]);
fver := fver * denv^2;
printf "FINALCURVE %o : y^2 = %o\n", Lab, fver;
Dfin := Z!Discriminant(fver);
// re-verify L-match on the final model, count split/inert (5|tp^2-4np pattern: inert <=> a_p in Z <=> tp^2-4np square)
nmatch := 0; nsplit := 0; ninert := 0;
for p in plist do
    if Dfin mod p eq 0 then continue; end if;
    tp, np := Explode(ap[p]);
    Ep := LpolyC(fver, p);
    if Ep ne TargetLpoly(p, tp, np) then
        printf "LMISMATCH_FINAL p=%o\n", p;
        continue;
    end if;
    nmatch +:= 1;
    if IsSquare(tp^2 - 4*np) then ninert +:= 1; else nsplit +:= 1; end if;
end for;
printf "LMATCH_FINAL %o : %o good primes match (deg2-split %o, deg1-rational %o)\n", Lab, nmatch, nsplit, ninert;

// ---- exact torsion ----
Jfin := Jacobian(HyperellipticCurve(fver));
t0 := Cputime();
T, mp := TorsionSubgroup(Jfin);
printf "TORSION_INVARIANTS %o : %o (%o s)\n", Lab, Invariants(T), Cputime(t0);
for i in [1..Ngens(T)] do
    P := mp(T.i);
    printf "TORSION_GEN %o: %o  order %o\n", i, P, Order(P);
end for;

// ---- simplicity certificate (D4 / root-power) ----
function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
    if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if;
    return cnt;
end function;
function SimplicityCertificate(fI)
    RT := PolynomialRing(Q); T := RT.1;
    dsc := Discriminant(fI);
    for pp in PrimesInInterval(3, 200) do
        if (Z!LeadingCoefficient(fI)) mod pp eq 0 then continue; end if;
        if (Z!Numerator(dsc)) mod pp eq 0 then continue; end if;
        PF := PolynomialRing(GF(pp));
        fp := PF![GF(pp)!co : co in Coefficients(fI)];
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        PF2 := PolynomialRing(GF(pp^2));
        fp2 := PF2![GF(pp^2)!co : co in Coefficients(fI)];
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
ok, pw, chiw := SimplicityCertificate(fver);
printf "SIMPLICITY %o : %o  witness prime %o  chi = %o\n", Lab, ok, pw, chiw;

// ---- repaired RM screen ----
cores := {};
for p in PrimesInInterval(3, 300) do
    if Dfin mod p eq 0 then continue; end if;
    Ep := LpolyC(fver, p);
    a1 := -Z!Coefficient(Ep, 1); a2 := Z!Coefficient(Ep, 2);
    n := a1^2 - 4*(a2 - 2*p);
    if a1 ne 0 and n ne 0 and not IsSquare(Z!Abs(n)) then Include(~cores, Squarefree(Z!Abs(n))); end if;
end for;
printf "RMSCREEN %o : distinct cores = %o => %o\n", Lab, Sort(Setseq(cores)), (#cores eq 1) select "RM" else "not-RM-shaped";
printf "BADPRIMES %o : %o (level primes %o)\n", Lab, BadPrimes(HyperellipticCurve(fver)), Factorization(N);
printf "VERIFYDONE %o\n", Lab;
quit;
