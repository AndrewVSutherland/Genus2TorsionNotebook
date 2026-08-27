// split_lab.m — shared infrastructure for the split-Jacobian torsion-table
// construction lanes (2026-08-12).  Load after SetColumns/SetMemoryLimit.
//
// KNOWN = groups already known to occur as J(Q)_tors of a geometrically split
// genus-2 Jacobian /Q with an *exact* witness.  Two-phase baseline:
//   - the 66 nontrivial groups (+ trivial) in the LMFDB alpha DB (g2c_curves_new,
//     is_simple_geom = false), queried 2026-08-12, plus
//   - repo exact split anchors: [45],[63],[70],[7,7] (results/claude_ov_lane8_verify.log),
//     [2,2,4,4] (notes/target_2248_source51_audit.md);
//   - PLUS, by default, the four groups realized by the 2026-08-12 session lanes
//     ([8,8],[6,12],[2,2,24],[2,2,4,8]; witnesses in data/new_split_witnesses.txt),
//     appended below so that FUTURE lane runs treat them as known.
// To reproduce the committed discovery logs (logs/glue2_*.log, x1fam.log,
// hlp37.log), which were produced BEFORE those four groups were known, set the
// environment variable SPLIT_KNOWN_BASELINE=2026-08-12-pre; the funnel then
// fires exact torsion on them again and emits the same HIT lines.

// Paths are checkout-relative: run every script in this directory FROM
// product/code/ (as the per-file Usage lines indicate).  The external
// dependency of the core lanes is Drew's Magma library (utils.m, polredabs.m,
// via its spec file), which is not part of this repository: set the
// environment variable TORSION_MAGMA_SPEC to its spec path, or rely on the
// default below.  (lane_z3.m additionally needs Frengley's N-congruences
// package; see that file's header.)
spec := GetEnv("TORSION_MAGMA_SPEC");
if spec eq "" then spec := "/home/claude/Magma/magma.spec"; end if;
AttachSpec(spec);
Attach("gluing.m");

KNOWN := [
  [Integers()|],
  [2],[3],[4],[5],[6],[7],[8],[9],[10],[12],[14],[15],[16],[18],[19],[20],[21],
  [24],[25],[27],[28],[30],[35],[36],[40],[48],[60],
  [2,2],[2,4],[2,6],[2,8],[2,10],[2,12],[2,14],[2,16],[2,18],[2,20],[2,24],[2,30],[2,48],
  [3,3],[3,6],[3,9],[3,12],[3,24],[4,4],[4,8],[4,12],[4,16],[5,5],[5,10],[6,6],
  [2,2,2],[2,2,4],[2,2,6],[2,2,8],[2,2,10],[2,2,12],[2,2,16],[2,4,4],[2,4,8],[2,6,6],
  [2,2,2,2],[2,2,2,4],[2,2,2,6],[2,2,2,8],
  [45],[63],[70],[7,7],[2,2,4,4]
];
if GetEnv("SPLIT_KNOWN_BASELINE") ne "2026-08-12-pre" then
    KNOWN cat:= [ [8,8],[6,12],[2,2,24],[2,2,4,8] ];   // realized 2026-08-12; see header
    KNOWN cat:= [ [11] ];   // realized 2026-08-26 (Q-simple conjugate-glue over Q(sqrt 11); new_split_witnesses.txt)
end if;

RQx<x> := PolynomialRing(Rationals());

function IntegralSextic(C)
    // y^2 = g(x) integral model of a genus-2 CrvHyp
    f,h := HyperellipticPolynomials(C);
    g := RQx!(4*f + h^2);
    if Degree(g) lt 5 or Degree(g) gt 6 then return false, g; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    g := RQx!(g*den^2);
    return true, g;
end function;

function OddPartInvs(seq)
    // canonical invariant chain of the odd part of an abelian group given by
    // any list of cyclic orders
    os := [ Integers() | d div 2^Valuation(d,2) : d in seq ];
    os := [ o : o in os | o ne 1 ];
    if IsEmpty(os) then return [Integers()|]; end if;
    return Invariants(AbelianGroup(os));
end function;

function ProfileCurve(g, NP)
    // invariants of J(F_p) at NP good odd primes (J(Q)_tors injects into each)
    dz := Integers()!Discriminant(g);
    lc := Integers()!LeadingCoefficient(g);
    bad := 2*dz*lc;
    invsList := []; ps := [];
    p := 2;
    while #invsList lt NP and p lt 1200 do
        p := NextPrime(p);
        if bad mod p eq 0 then continue; end if;
        gp := PolynomialRing(GF(p))!g;
        Jp := Jacobian(HyperellipticCurve(gp));
        Ap := AbelianGroup(Jp);
        Append(~invsList, Invariants(Ap)); Append(~ps, p);
    end while;
    return invsList, ps;
end function;

function EllProfile(invsList, ell)
    // a_k = min over p of #G_p[ell^k] as an exponent, k = 1,2,... until stable
    vexps := [ [ Valuation(d, ell) : d in invs ] : invs in invsList ];
    aks := [Integers()|];
    kk := 1;
    repeat
        ak := Min([ IsEmpty(ve) select 0 else &+[ Min(kk, v) : v in ve ] : ve in vexps ]);
        Append(~aks, ak);
        stable := kk gt 1 and aks[kk] eq aks[kk-1];
        kk +:= 1;
    until stable or kk gt 14;
    return aks;
end function;

function EllShapes(aks, ell)
    // all exponent shapes m1>=m2>=m3>=m4 (<=2 nonzero parts if ell odd, since
    // rational odd ell-torsion of an abelian surface has rank <= 2) with
    // sum_i min(m_i,k) <= a_k for all k
    maxr := ell eq 2 select 4 else 2;
    M := aks[#aks];
    shapes := [];
    for m1 in [0..M] do
    for m2 in [0..(maxr ge 2 select m1 else 0)] do
    for m3 in [0..(maxr ge 3 select m2 else 0)] do
    for m4 in [0..(maxr ge 4 select m3 else 0)] do
        ms := [m1,m2,m3,m4];
        if &+ms gt M then continue; end if;
        ok := true;
        for k in [1..#aks] do
            if &+[ Min(ms[i],k) : i in [1..4] ] gt aks[k] then ok := false; break; end if;
        end for;
        if ok then Append(~shapes, ms); end if;
    end for; end for; end for; end for;
    return shapes;
end function;

function CompatibleGroups(invsList)
    // all abelian groups T that embed into every J(F_p) profile
    ords := [ IsEmpty(invs) select 1 else &*invs : invs in invsList ];
    d := GCD(ords);
    if d eq 1 then return [ [Integers()|] ]; end if;
    combos := [ [Parent(<2,[1,0,0,0]>)|] ];
    for ell in PrimeDivisors(d) do
        aks := EllProfile(invsList, ell);
        shs := EllShapes(aks, ell);
        combos := [ c cat [<ell,s>] : c in combos, s in shs ];
    end for;
    out := [];
    for c in combos do
        chain := [1,1,1,1];
        for t in c do
            for j in [1..4] do chain[j] *:= t[1]^(t[2][j]); end for;
        end for;
        inv := Reverse([ x : x in chain | x ne 1 ]);
        Append(~out, inv);
    end for;
    return Setseq(Seqset(out));
end function;

function ExactTorsion(g)
    // TorsionSubgroup requires a y^2 = f(x) integral model
    okm := false; gm := g;
    try
        Cm := ReducedMinimalWeierstrassModel(HyperellipticCurve(g));
        okm, gm := IntegralSextic(Cm);
    catch e
        okm := false;
    end try;
    if okm then g := gm; end if;
    J := Jacobian(HyperellipticCurve(g));
    return Invariants(TorsionSubgroup(J));
end function;

function LpReducibleCount(g, pmax)
    dz := Integers()!Discriminant(g);
    lc := Integers()!LeadingCoefficient(g);
    bad := 2*dz*lc;
    n := 0; red := 0;
    p := 2;
    while p lt pmax do
        p := NextPrime(p);
        if bad mod p eq 0 then continue; end if;
        gp := PolynomialRing(GF(p))!g;
        Lp := Numerator(ZetaFunction(HyperellipticCurve(gp)));
        n +:= 1;
        if not IsIrreducible(Lp) then red +:= 1; end if;
    end while;
    return red, n;
end function;

function Funnel(C, tag : OddInvs := [Integers()|-1])
    // returns one of: "skip","abort","known","hit","fail".
    // OddInvs: if not [-1], the exactly-known odd part of J(Q)_tors (e.g. for
    // (2,2)-glued curves it is T1_odd x T2_odd); used to prune candidates.
    ok, g := IntegralSextic(C);
    if not ok then printf "SKIP %o deg\n", tag; return "skip"; end if;
    if Discriminant(g) eq 0 then printf "SKIP %o sing\n", tag; return "skip"; end if;
    invsList, ps := ProfileCurve(g, 12);
    if #invsList lt 4 then printf "SKIP %o fewprimes\n", tag; return "skip"; end if;
    compat := CompatibleGroups(invsList);
    if OddInvs ne [Integers()|-1] then
        compat := [ v : v in compat | OddPartInvs(v) eq OddInvs ];
    end if;
    unknown := [ v : v in compat | not v in KNOWN ];
    if #unknown eq 0 then
        return "abort";
    end if;
    printf "CAND %o unknown=%o\n", tag, unknown;
    einv := [Integers()|-1];
    try
        einv := ExactTorsion(g);
    catch e
        printf "EXACTFAIL %o g=%o\n", tag, g;
        return "fail";
    end try;
    if einv in KNOWN then
        printf "EXACTKNOWN %o invs=%o\n", tag, einv;
        return "known";
    end if;
    red, n := LpReducibleCount(g, 120);
    printf "HIT %o invs=%o Lpred=%o/%o g=%o\n", tag, einv, red, n, g;
    return "hit";
end function;
