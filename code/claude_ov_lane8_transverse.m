// Lane 8 (2026-07-25, resumed session): TRANSVERSE-DEFORMATION experiment.
//
// Lane item 1: "put the split point in a chart, write the torsion condition
// around it, identify the split divisor through it, and determine whether the
// torsion locus has a component TRANSVERSE to the split locus."
//
// Geometry.  Every recovered HLP model is bielliptic: y^2 = a x^6+b x^4+c x^2+d,
// with the extra involution x -> -x.  The bielliptic (= (2,2)-split) locus L in
// the 3-dimensional M_2 is a SURFACE: 4 even coefficients minus the 1-parameter
// reparametrization x -> lambda x minus the y-scaling.  So codim L = 1 and the
// transverse direction is 1-dimensional (the 3 odd coefficients modulo the 2
// PGL_2 directions that do not preserve the involution).
//
// Note on infinitesimal deformation: in characteristic 0 the group scheme J[N]
// is FINITE ETALE over the base, so a rational N-torsion section lifts UNIQUELY
// over any Artinian thickening.  First-order deformation theory therefore gives
// NO obstruction and NO information.  The only meaningful test is a global one:
// does the torsion condition hold along a positive-dimensional family through
// the anchor?  If it did over Q, it would hold for EVERY good t in F_p.
//
// So: reduce mod p, walk the transverse lines through the anchor, and measure
// the density of the target-torsion condition on chi_t against the density over
// random curves at the same p (the "baseline").
SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned PMAX then PMAX := 400; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;
if not assigned NRAND then NRAND := 4; elif Type(NRAND) eq MonStgElt then NRAND := StringToInteger(NRAND); end if;
SetMemoryLimit(MemGB*10^9);
QQ := Rationals();
P<x> := PolynomialRing(QQ);

// --- the recovered HLP anchors (from results/claude_ov_lane8_verify.log) ---
anchors := [
  <"Z/5xZ/10", 45600*x^6 + 289161*x^4 + 35186670*x^2 - 705688215, 50, [<5,2>], 2>,
  <"Z/7xZ/7",  x^6 + 3025*x^4 + 3232987*x^2 + 869675859,          49, [<7,2>], 1>,
  <"Z/63",     897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496,   63, [],      1>,
  <"Z/45",     13981*x^6 + 29240200*x^4 + 49996210000*x^2 + 168300000000, 45, [], 1>,
  <"Z/35",     640*x^6 + 5040*x^4 + 2480*x^2 + 9295,              35, [],      1>
];
// entry: <tag, f0 (even sextic), N = |G|, list of <n,2> rank-2 conditions, 2-rank needed>

// condition on the Frobenius charpoly chi of a genus-2 curve over F_p
function CompatOK(chi, p, N, ranks)
    nJ := Integers()!Evaluate(chi, 1);
    if nJ mod N ne 0 then return false; end if;
    for r in ranks do
        n := r[1];
        if nJ mod (n^2) ne 0 then return false; end if;
        Rn := PolynomialRing(GF(n));
        if (Rn!chi) mod (Rn.1-1)^2 ne 0 then return false; end if;
    end for;
    return true;
end function;

function ChiOf(fp)
    if Degree(fp) lt 5 or not IsSquarefree(fp) then return 0, false; end if;
    ok := true; C := 0;
    try C := HyperellipticCurve(fp); catch e ok := false; end try;
    if not ok then return 0, false; end if;
    if Genus(C) ne 2 then return 0, false; end if;
    L := LPolynomial(C);
    // LPolynomial is the reciprocal: L(T) = prod(1-alpha T).  chi(T) = T^4 L(1/T)
    RZ := PolynomialRing(Integers()); T := RZ.1;
    cs := Coefficients(L);
    chi := &+[ cs[i]*T^(5-i) : i in [1..#cs] ];
    return chi, true;
end function;

for A in anchors do
    tag := A[1]; f0 := A[2]; N := A[3]; ranks := A[4];
    D0 := Integers()!Discriminant(f0);
    bad := {q[1] : q in Factorization(D0)} join {2};
    printf "\n############ %o  anchor y^2 = %o\n", tag, f0;
    printf "  bad primes of the anchor : %o\n", Sort(Setseq(bad));
    totline := 0; hitline := 0; totrand := 0; hitrand := 0; totflat := 0; hitflat := 0;
    allgood := true;
    for p in PrimesUpTo(PMAX) do
        if p in bad then continue; end if;
        Fp := GF(p); Pp := PolynomialRing(Fp); X := Pp.1;
        fp0 := Pp!f0;
        // sanity: does the anchor itself pass at p?
        chi0, ok0 := ChiOf(fp0);
        if not ok0 or not CompatOK(chi0, p, N, ranks) then allgood := false; end if;
        // (1) the three transverse lines  f0 + t*x^k, k = 1,3,5
        for k in [1,3,5] do
            for t in Fp do
                if t eq 0 then continue; end if;
                ft := fp0 + t*X^k;
                chi, ok := ChiOf(ft);
                if not ok then continue; end if;
                totline +:= 1;
                if CompatOK(chi, p, N, ranks) then hitline +:= 1; end if;
            end for;
        end for;
        // (2) the two "flat" (bielliptic-preserving) lines f0 + t*x^4, t*x^2 : stay ON L
        for k in [0,2,4] do
            for t in Fp do
                if t eq 0 then continue; end if;
                ft := fp0 + t*X^k;
                chi, ok := ChiOf(ft);
                if not ok then continue; end if;
                totflat +:= 1;
                if CompatOK(chi, p, N, ranks) then hitflat +:= 1; end if;
            end for;
        end for;
        // (3) baseline: random sextics over F_p
        for i in [1..NRAND*p] do
            ft := Pp![Random(Fp) : j in [1..7]];
            if Degree(ft) lt 5 then continue; end if;
            chi, ok := ChiOf(ft);
            if not ok then continue; end if;
            totrand +:= 1;
            if CompatOK(chi, p, N, ranks) then hitrand +:= 1; end if;
        end for;
    end for;
    printf "  anchor passes the mod-p compatibility test at every good p < %o : %o\n", PMAX, allgood;
    printf "  TRANSVERSE lines  f0 + t*x^{1,3,5} : %o / %o = %o\n",
           hitline, totline, totline eq 0 select 0 else RealField(5)!(hitline/totline);
    printf "  IN-LOCUS   lines  f0 + t*x^{0,2,4} : %o / %o = %o\n",
           hitflat, totflat, totflat eq 0 select 0 else RealField(5)!(hitflat/totflat);
    printf "  BASELINE   random sextics          : %o / %o = %o\n",
           hitrand, totrand, totrand eq 0 select 0 else RealField(5)!(hitrand/totrand);
end for;

printf "SEARCH_DONE transverse\n";
quit;
