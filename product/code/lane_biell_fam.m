// lane_biell_fam.m — bielliptic x -> x^2 sweep over the GENUS-0 TORSION
// FAMILIES (extends lane_biell.m, which only used LMFDB seeds; AVS request
// 2026-08-26): for E1 = E1(t) in a Kubert/HLP family with interesting
// torsion, put E1 in the form y^2 = f_t(x), shift g(x) := f_t(x+s), and
// funnel the genus-2 curve
//     C_{t,s} : y^2 = g(x^2),
// whose Jacobian is the (2,2)-gluing of E1 (u = x^2) with the "reversed"
// partner E2 : v^2 = g0 u^3 + g1 u^2 + g2 u + g3 (u = 1/x^2, v = y/x^3).
// The shift s moves E2 over a 1-parameter family of 2-congruent partners of
// E1.  Odd torsion of the glue is exactly T1_odd x T2_odd; 2-primary parts
// can gain.  Funnel gates:
//   (i) E2 has torsion of order >= 5 (checked exactly after a mod-p gcd
//       prefilter), or
//  (ii) a 2-primary gain is possible: v2(#E1(F_p) * #E2(F_p)) >=
//       v2(|T1[2^inf]| * |T2[2^inf]|) + 1 at ALL test primes.
// Survivors go to split_lab's Funnel (early abort vs KNOWN before any exact
// genus-2 torsion computation), with OddInvs pinned to T1_odd x T2_odd.
//
// Usage (from product/code/):
//   magma -b Fam:=12 TH:=24 SH:=24 Part:=1 NParts:=1 lane_biell_fam.m > ../logs/biellfam_12_p1.log
// Fam in {5,6,7,8,9,10,12,24,26,28} (24=[2,4], 26=[2,6], 28=[2,8]).
SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned Fam then Fam := 12; elif Type(Fam) eq MonStgElt then Fam := StringToInteger(Fam); end if;
if not assigned TH then TH := 24; elif Type(TH) eq MonStgElt then TH := StringToInteger(TH); end if;
if not assigned SH then SH := 24; elif Type(SH) eq MonStgElt then SH := StringToInteger(SH); end if;
if not assigned Part then Part := 1; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;

load "split_lab.m";  // run from product/code/

Q := Rationals();
Px<x> := PolynomialRing(Q);

// families as monic-cleared cubics y^2 = f(x); Tate form completed to
// (2y+a1x+a3)^2 = 4x^3 + (a1^2+4a2)x^2 + 2a1a3 x + a3^2
function TateCubic(bv, cv)
    // build the Tate-form curve directly and return ok + the MONIC cubic of
    // its short Weierstrass model (EllipticCurve requires monic cubics)
    E := 0;
    try E := EllipticCurve([1-cv, -bv, -bv, 0, 0]); catch e; return false, x; end try;
    S := WeierstrassModel(E);
    a := aInvariants(S);
    return true, x^3 + a[4]*x + a[5];
end function;

function FamilyCubic(F, tv)
    // returns ok, f (cubic in x with y^2 = f defining the family curve)
    if F eq 5 then
        return TateCubic(tv, tv);
    elif F eq 6 then
        return TateCubic(tv^2+tv, tv);
    elif F eq 7 then
        return TateCubic(tv^3-tv^2, tv^2-tv);
    elif F eq 8 then
        if tv eq 0 then return false, x; end if;
        return TateCubic((2*tv-1)*(tv-1), (2*tv-1)*(tv-1)/tv);
    elif F eq 9 then
        return TateCubic(tv^2*(tv-1)*(tv^2-tv+1), tv^2*(tv-1));
    elif F eq 10 then
        de := tv^2-3*tv+1;
        if de eq 0 then return false, x; end if;
        return TateCubic(tv^3*(tv-1)*(2*tv-1)/de^2, -tv*(tv-1)*(2*tv-1)/de);
    elif F eq 12 then
        if tv eq 1 then return false, x; end if;
        return TateCubic(tv*(2*tv-1)*(3*tv^2-3*tv+1)*(2*tv^2-2*tv+1)/(tv-1)^4,
                               -tv*(2*tv-1)*(3*tv^2-3*tv+1)/(tv-1)^3);
    elif F eq 24 then
        return true, x*(x+1)*(x+tv^2);
    elif F eq 26 then
        if (tv+3)*(tv-3) eq 0 then return false, x; end if;
        r1 := (-2*tv+10)/((tv+3)*(tv-3));
        r2 := (-tv^3+7*tv^2-11*tv+5)/(4*(tv+3)*(tv-3)^2);
        r3 := (-2*tv^2+4*tv-2)/((tv+3)^2*(tv-3));
        return true, (x-r1)*(x-r2)*(x-r3);
    elif F eq 28 then
        if tv eq 0 or 8*tv^2-1 eq 0 or 4*tv+1 eq 0 then return false, x; end if;
        r1 := (16*tv^3+12*tv^2+2*tv)/(8*tv^2-1)^2;
        r2 := (32*tv^3+24*tv^2+8*tv+1)/(16*tv^2*(8*tv^2-1));
        r3 := (-32*tv^4-32*tv^3-12*tv^2-2*tv)/((4*tv+1)^2*(8*tv^2-1));
        return true, (x-r1)*(x-r2)*(x-r3);
    else
        error "unknown family", F;
    end if;
end function;

EXPECTED := AssociativeArray();
EXPECTED[5] := [5]; EXPECTED[6] := [6]; EXPECTED[7] := [7]; EXPECTED[8] := [8];
EXPECTED[9] := [9]; EXPECTED[10] := [10]; EXPECTED[12] := [12];
EXPECTED[24] := [2,4]; EXPECTED[26] := [2,6]; EXPECTED[28] := [2,8];

// runtime family verification
nv := 0;
for tv in [Q| 3, 5/2, -7/3 ] do
    ok, f := FamilyCubic(Fam, tv);
    if not ok or Discriminant(f) eq 0 then continue; end if;
    try
        if Invariants(TorsionSubgroup(EllipticCurve(f))) eq EXPECTED[Fam] then nv +:= 1; end if;
    catch e;
    end try;
end for;
error if nv lt 2, "family", Fam, "failed runtime verification";
printf "FAMILY %o verified %o/3\n", Fam, nv;

T1inv := EXPECTED[Fam];
T1odd := [ i div 2^Valuation(i,2) : i in T1inv | i div 2^Valuation(i,2) gt 1 ];
v2T1 := &+[ Valuation(i, 2) : i in T1inv ];

PRE_P := [101, 103, 107, 109, 113, 127, 131, 137];

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, Q!a/b); Append(~S, Q!-a/b); end if;
    end for;
    return S;
end function;

function MonicModel(cub)
    // y^2 = a x^3 + b x^2 + c x + d  ~  Y^2 = X^3 + b X^2 + ac X + a^2 d
    a := Coefficient(cub,3); b := Coefficient(cub,2);
    c := Coefficient(cub,1); d := Coefficient(cub,0);
    return x^3 + b*x^2 + a*c*x + a^2*d;
end function;

function CountCubic(f, p)
    // #E(F_p) for y^2 = f (cubic); false when f doesn't reduce well mod p
    Fp := GF(p);
    for c in Coefficients(f) do
        if Denominator(c) mod p eq 0 then return 0, false; end if;
    end for;
    Rp := PolynomialRing(Fp);
    fp := Rp![ Fp!Numerator(c)/Fp!Denominator(c) : c in Coefficients(f) ];
    if Degree(fp) ne 3 or Discriminant(fp) eq 0 then return 0, false; end if;
    return #EllipticCurve(fp), true;
end function;

tvals := HeightRats(TH);
svals := [Q|0] cat HeightRats(SH);
nt := 0; npair := 0; ngate1 := 0; ngate2 := 0; nfun := 0; nhit := 0; nabort := 0;
t0c := Cputime();

for ti in [1..#tvals] do
    if ti mod NParts ne (Part - 1) then continue; end if;
    tv := tvals[ti];
    ok, f := FamilyCubic(Fam, tv);
    if not ok or Discriminant(f) eq 0 then continue; end if;
    okT := false;
    try okT := Invariants(TorsionSubgroup(EllipticCurve(f))) eq T1inv; catch e; end try;
    if not okT then continue; end if;   // stay on the exact stratum
    nt +:= 1;
    // cache E1 counts
    n1 := AssociativeArray();
    for p in PRE_P do
        c, okc := CountCubic(f, p);
        if okc then n1[p] := c; end if;
    end for;
    for sv in svals do
        g := Evaluate(f, x + sv);
        g0 := Coefficient(g, 0);
        if g0 eq 0 then continue; end if;
        npair +:= 1;
        rev := g0*x^3 + Coefficient(g,1)*x^2 + Coefficient(g,2)*x + Coefficient(g,3);
        revm := MonicModel(rev);
        // E2 counts at the test primes (cached for both gates)
        c2s := AssociativeArray();
        for p in PRE_P do
            if not IsDefined(n1, p) then continue; end if;
            c2, okc := CountCubic(revm, p);
            if okc then c2s[p] := c2; end if;
        end for;
        ks := [ p : p in PRE_P | IsDefined(c2s, p) ];
        if #ks lt 5 then continue; end if;
        gcdv := 0;
        for p in ks do gcdv := GCD(gcdv, c2s[p]); end for;
        // ---- gate (i): E2 torsion order >= 5?  (exact check whenever the
        // gcd leaves room for order >= 5 OR for odd 3-torsion, which the
        // OddInvs pin below must know about)
        T2inv := [Integers()|];
        if gcdv ge 5 or gcdv eq 3 then
            try T2inv := Invariants(TorsionSubgroup(EllipticCurve(revm))); catch e; end try;
        end if;
        ord2 := #T2inv eq 0 select 1 else &*T2inv;
        gate1 := ord2 ge 5;
        // ---- gate (ii): 2-primary gain possible?
        gate2 := true;
        if not gate1 then
            if gcdv ge 5 or gcdv eq 3 then
                v2T2 := #T2inv eq 0 select 0 else &+[ Valuation(i,2) : i in T2inv ];
            else
                nr := #Roots(rev);
                v2T2 := nr ge 1 select (nr eq 3 select 2 else 1) else 0;
            end if;
            need := v2T1 + v2T2 + 1;
            for p in ks do
                if Valuation(n1[p]*c2s[p], 2) lt need then gate2 := false; break; end if;
            end for;
        end if;
        if not (gate1 or gate2) then continue; end if;
        if gate1 then ngate1 +:= 1; else ngate2 +:= 1; end if;
        // ---- funnel the genus-2 curve
        C6 := 0; okC := true;
        try
            sext := Evaluate(g, x^2);
            if Discriminant(sext) eq 0 then okC := false; else C6 := HyperellipticCurve(sext); end if;
        catch e;
            okC := false;
        end try;
        if not okC then continue; end if;
        T2odd := [ i div 2^Valuation(i,2) : i in T2inv | i div 2^Valuation(i,2) gt 1 ];
        oddprod := Sort([ Integers() | i : i in T1odd cat T2odd ]);
        // normalize odd product to invariant factors
        oddinvs := #oddprod eq 0 select [Integers()|] else Invariants(AbelianGroup(oddprod));
        nfun +:= 1;
        res := Funnel(C6, Sprintf("biellfam|%o(t=%o)s=%o|T2=%o", Fam, tv, sv, T2inv)
                      : OddInvs := oddinvs);
        if res eq "hit" then nhit +:= 1; end if;
        if res eq "abort" then nabort +:= 1; end if;
    end for;
    if nt mod 25 eq 0 then
        printf "PROGRESS t#%o pairs=%o g1=%o g2=%o fun=%o abort=%o hit=%o %os\n",
            nt, npair, ngate1, ngate2, nfun, nabort, nhit, Cputime()-t0c;
    end if;
end for;
printf "SEARCH_DONE biellfam Fam=%o TH=%o SH=%o part %o/%o t=%o pairs=%o gate1=%o gate2=%o funneled=%o aborts=%o hits=%o %.1o s\n",
    Fam, TH, SH, Part, NParts, nt, npair, ngate1, ngate2, nfun, nabort, nhit, Cputime()-t0c;
quit;
