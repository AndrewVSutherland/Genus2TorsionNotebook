//////////////////////////////////////////////////////////////////////
// The elliptic quotient of the Prym over K = Q(zeta_3).
// Prym ~ J(C2'), C2'_min: y^2 = -3x^6+24x^3-75; its (-3)-twist is
// D: y^2 = x^6 - 8x^3 + 25  (nicer; has rational points).
// mu_3 structure => elliptic quotients over K have j = 0:
//   E_d: v^2 = u^3 + d,  d in K* / (K*)^6,  J(D)_K ~ E_d x E_d^sigma.
// This script: (1) rank of J(D) (the twist rank r_tw);
// (2) exact quadratic splitting of L_p(J(D)) at split primes;
// (3) scan sextic-twist candidates d, match a_p pairs, output d.
// Usage: magma -b code/contact6_m612_relative3_prym_elliptic.m
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1); SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);
fD := x^6 - 8*x^3 + 25;
D := HyperellipticCurve(fD);
JD := Jacobian(D);
rlo, rhi := RankBounds(JD);
printf "TWIST RANK BOUNDS (J of y^2=x^6-8x^3+25): %o..%o\n", rlo, rhi;
T := TorsionSubgroup(JD);
printf "twist torsion: %o\n", Invariants(T);

// L-factors at split primes p = 1 mod 3
splitp := [7,13,19,31,37,43,61,67];
PT<Tt> := PolynomialRing(Integers());
Lfacs := AssociativeArray();
for p in splitp do
    Fp := GF(p);
    fp := PolynomialRing(Fp)!fD;
    if not IsSquarefree(fp) then continue; end if;
    Dp := HyperellipticCurve(fp);
    Zf := Numerator(ZetaFunction(Dp));   // degree 4, = L_p(T)
    fac := Factorization(PT!Zf);
    printf "p=%o : L_p factors %o\n", p, fac;
    Lfacs[p] := fac;
end for;

// candidate scan for d: units {+-1, +-z, +-z^2} * 2^a * lam^b * 5^c, lam = zeta-zeta^2 = sqrt(-3)
K<z> := CyclotomicField(3);
lam := z - z^2;      // sqrt(-3)
OK := Integers(K);
units := [K | 1, -1, z, -z, z^2, -z^2];
cands := [];
for u in units do for a in [0..5] do for b in [0..5] do for c in [0..5] do
    Append(~cands, u * 2^a * lam^b * 5^c);
end for; end for; end for; end for;
printf "candidates: %o\n", #cands;

// for each split p: the two embeddings K -> F_p; a_p of E_d mod each
function APair(d, p)
    Fp := GF(p);
    rts := [r[1] : r in Roots(PolynomialRing(Fp)![1,1,1])];  // z -> rts
    assert #rts eq 2;
    aps := [];
    for r in rts do
        // evaluate d at z -> r
        dp := &+[Fp | Fp!(Integers()!Eltseq(d)[i]) * r^(i-1) : i in [1..2]]
              where _ := 0;
        if dp eq 0 then return false, []; end if;
        Ep := EllipticCurve([Fp | 0, dp]);
        Append(~aps, p + 1 - #Ep);
    end for;
    return true, Sort(aps);
end function;

// required a-pairs from the L-factor quadratic splitting:
// L_p = (1 - a T + p T^2)(1 - a' T + p T^2): trace coefficients
function ReqPairs(p, Lfac)
    // each degree-2 factor T^2 - a T + p ... (as poly in T with Magma's
    // zeta numerator convention: factors are p*T^2 - a*T + 1 up to ordering)
    prs := [];
    if #Lfac eq 2 and Degree(Lfac[1][1]) eq 2 then
        a1 := -Coefficient(Lfac[1][1], 1) / Coefficient(Lfac[1][1], 0) * 1;
        a2 := -Coefficient(Lfac[2][1], 1) / Coefficient(Lfac[2][1], 0) * 1;
        // normalize: constant coeff 1 convention: a = -c1/c0
        Append(~prs, Sort([Integers()!a1, Integers()!a2]));
    end if;
    return prs;
end function;

survivors := cands;
for p in splitp do
    if not IsDefined(Lfacs, p) then continue; end if;
    fac := Lfacs[p];
    if #fac ne 2 then printf "p=%o: not split into 2 factors; skip\n", p; continue; end if;
    req := ReqPairs(p, fac);
    if #req eq 0 then continue; end if;
    newsurv := [];
    for d in survivors do
        ok, ap := APair(d, p);
        if ok and ap eq req[1] then Append(~newsurv, d); end if;
    end for;
    survivors := newsurv;
    printf "p=%o : required a-pair %o ; survivors now %o\n", p, req[1], #survivors;
    if #survivors le 12 then break; end if;
end for;
printf "SURVIVING d candidates: %o\n", survivors;
print "PRYM_ELLIPTIC_DONE";
quit;
