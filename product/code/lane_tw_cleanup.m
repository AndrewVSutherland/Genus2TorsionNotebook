// lane_tw_cleanup.m — Stage 5b: make the twisted-diagonal impossibility
// unconditional (session-note 2026-08-13 item 2.4 -> 2026-08-14 plan 3.3).
// The two genus-3 twisted condition curves (exactly lane_misc2.m:19-20):
//   tw12: w^2 = -1296*T^8 + 5184*T^7 - 9072*T^6 + 9072*T^5 - 5580*T^4
//               + 2088*T^3 - 432*T^2 + 36*T
//   tw10: w^2 = 32*T^7 - 160*T^6 + 256*T^5 - 156*T^4 + 16*T^3 + 16*T^2 - 4*T
// lane_misc2.m established RankBound(J) = 0 for both, under
// SetClassGroupBounds("GRH"), and listed the small rational points.  This lane
// completes the argument:
//   (1) pass to an ODD-degree monic integral model (tw12: T -> 1/x flip, strip
//       the square content 36 = 6^2; both: monicize by (x,y) -> (cx, c^3 y));
//   (2) TORSION BOUND: #J(F_p) = L_p(1) at >= 6 good odd primes; gcd -> B.
//       (Reduction J(Q)_tors -> J(F_p) is injective for odd p of good
//       reduction, so B bounds #J(Q)_tors; several primes make gcd tight.)
//   (3) KNOWN SUBGROUP: classes [P - infty] for every small rational point P
//       (rational Weierstrass roots + Points(C : Bound := 10^4)); close under
//       addition -> subgroup S <= J(Q).
//   (4) If #S = B and rank 0: J(Q) = S exactly.  Enumerate S; on an odd model
//       reduced Mumford representatives are unique, and an element is [P-infty]
//       for P in C(Q) iff its a-polynomial has degree <= 1 (identity <-> infty).
//       Since C(Q) -> J(Q), P |-> [P - infty], is injective (g >= 1), C(Q) is
//       EXACTLY {infty} u {degree-1 elements}: prints CLEANUP ... CQ_EXHAUSTED.
//   (5) LAST (Alarm-guarded): attempt RankBound WITHOUT GRH class-group bounds
//       for a fully unconditional rank 0; SKIPPED/FAILED just means the rank-0
//       input stays GRH-conditional (say so in the writeup).
// NOTE on Alarm: Alarm(AlarmS) kills the process if the unconditional
// RankBound hangs; -b buffering could then lose the stdout tail, so every key
// result is ALSO mirrored immediately to ../logs/twcleanup.progress.
// Markers: TWCLEANUP_START, TORSBOUND, KNOWNSUB, CLEANUP, RANKBOUND_UNCOND,
//          TWCLEANUP_DONE
// Usage: cd product/code && magma -b lane_tw_cleanup.m > ../logs/lane_tw_cleanup.log
//   optional: MemGB:=<int> (default 5), AlarmS:=<int> (default 1800),
//             DoRank:=<0|1> (default 1)
SetColumns(0);
if not assigned MemGB then MemGB := 5; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned AlarmS then AlarmS := 1800; elif Type(AlarmS) eq MonStgElt then AlarmS := StringToInteger(AlarmS); end if;
if not assigned DoRank then DoRank := 1; elif Type(DoRank) eq MonStgElt then DoRank := StringToInteger(DoRank); end if;
SetMemoryLimit(MemGB*10^9);
SetSeed(1);
Alarm(AlarmS);

RQ := Rationals();
PR<x> := PolynomialRing(RQ);

function Monicize(f)   // odd-degree f -> monic integral F; (x,y) -> (c x, c^((d-1) div 2) y)
    d := Degree(f); c := LeadingCoefficient(f);
    F := PR ! &+[ Coefficient(f, i) * c^(d-1-i) * x^i : i in [0..d] ];
    return F, c;
end function;

function StripSquareContent(g)   // g = s^2 * h with h integral of trivial square content
    cont := GCD([ Integers() ! c : c in Coefficients(g) ]);
    _, s := SquarefreeFactorization(cont);
    return PR ! (g / s^2), s;
end function;

printf "TWCLEANUP_START AlarmS=%o DoRank=%o\n", AlarmS, DoRank;
System("rm -f ../logs/twcleanup.progress");

// ---- build the odd monic models ----
f12 := -1296*x^8 + 5184*x^7 - 9072*x^6 + 9072*x^5 - 5580*x^4 + 2088*x^3 - 432*x^2 + 36*x;
f10 := 32*x^7 - 160*x^6 + 256*x^5 - 156*x^4 + 16*x^3 + 16*x^2 - 4*x;

// tw12: T | f12 and deg 8, so T -> 1/x is an odd-degree model: g(x) = x^8 f12(1/x)
g12 := PR ! (x^8 * Evaluate(f12, 1/x));
h12, s12 := StripSquareContent(g12);           // expect s12 = 6, h12 monic deg 7
F12, c12 := Monicize(h12);
error if Degree(F12) ne 7 or LeadingCoefficient(F12) ne 1, "tw12 model construction failed";
// point map tw12: (T,w) with T ne 0  ->  x = c12/T (c12 = 1 expected), and T=0 -> infty
F10, c10 := Monicize(f10);                     // point map: (T,w) -> (32 T, 32^3 w); infty -> infty
error if Degree(F10) ne 7 or LeadingCoefficient(F10) ne 1, "tw10 model construction failed";
printf "MODEL tw12: y^2 = %o  (T = %o/x; infty <-> T=0)\n", F12, c12;
printf "MODEL tw10: y^2 = %o  (T = x/%o; infty <-> T=infty)\n", F10, c10;

procedure DoCurve(name, F, tlabel)
    error if Discriminant(F) eq 0, "model not squarefree";
    C := HyperellipticCurve(F);
    J := Jacobian(C);
    printf "== %o: genus %o ==\n", name, Genus(C);
    error if Genus(C) ne 3, "genus drift";

    // (2) torsion bound via #J(F_p) = L_p(1) at good odd primes
    D := Integers() ! Discriminant(F);
    ps := []; p := 3;
    while #ps lt 6 do
        if D mod p ne 0 then Append(~ps, p); end if;
        p := NextPrime(p);
    end while;
    B := 0;
    for p in ps do
        n := 0;
        try
            Cp := HyperellipticCurve(PolynomialRing(GF(p)) ! F);
            n := Integers() ! Evaluate(Numerator(ZetaFunction(Cp)), 1);
        catch e
            n := 0;
        end try;
        if n eq 0 then
            try n := #BaseChange(J, GF(p)); catch e n := 0; end try;
        end if;
        if n ne 0 then
            printf "  #J(F_%o) = %o = %o\n", p, n, Factorization(n);
            B := Gcd(B, n);
        else
            printf "  #J(F_%o) FAILED\n", p;
        end if;
    end for;
    printf "TORSBOUND curve=%o B=%o primes=%o\n", name, B, ps;
    System(Sprintf("echo 'TORSBOUND %o B=%o' >> ../logs/twcleanup.progress", name, B));

    // (3) known rational points -> classes [P - infty] -> subgroup closure
    knownpts := [];   // affine [x0, y0]
    for r in Roots(F) do Append(~knownpts, [r[1], RQ!0]); end for;
    try
        for P in Points(C : Bound := 10^4) do
            if P[3] ne 0 then
                a := [P[1]/P[3], P[2]/P[3]^4];
                if not a in knownpts then Append(~knownpts, a); end if;
            end if;
        end for;
    catch e
        printf "  point search failed (Weierstrass roots only)\n";
    end try;
    printf "  known points: infty + %o\n", knownpts;
    gens := [];
    for a in knownpts do
        g := J ! 0; okg := false;
        try g := elt<J | x - a[1], PR ! a[2], 1>; okg := true; catch e; end try;
        if not okg then
            try g := J ! [x - a[1], PR ! a[2]]; okg := true; catch e; end try;
        end if;
        if not okg then
            try
                g := C ! [a[1], a[2]] - PointsAtInfinity(C)[1];
                okg := true;
            catch e; end try;
        end if;
        error if not okg, Sprintf("class construction failed for %o on %o", a, name);
        Append(~gens, g);
    end for;
    S := {@ J ! 0 @};
    frontier := [ J ! 0 ];
    while #frontier gt 0 and #S le 512 do
        newelts := [];
        for s in frontier, g in gens do
            e := s + g;
            if not e in S then Include(~S, e); Append(~newelts, e); end if;
        end for;
        frontier := newelts;
    end while;
    printf "KNOWNSUB curve=%o order=%o\n", name, #S;
    System(Sprintf("echo 'KNOWNSUB %o order=%o' >> ../logs/twcleanup.progress", name, #S));

    // (4) exhaustiveness
    if #S eq B then
        npts := 0;
        for e in S do
            d := Degree(e[1]);
            if d eq 0 then
                printf "  element 0 <-> infty (%o)\n", tlabel[1];
                npts +:= 1;
            elif d eq 1 then
                x0 := -Coefficient(e[1], 0)/Coefficient(e[1], 1);
                y0 := Evaluate(e[2], x0);
                printf "  element deg1 <-> point (%o, %o), T = %o\n", x0, y0,
                    tlabel[2] eq "inv" select tlabel[3]/x0 else x0/tlabel[3];
                npts +:= 1;
            else
                printf "  element deg%o (not a point class)\n", d;
            end if;
        end for;
        printf "CLEANUP curve=%o CQ_EXHAUSTED npts=%o (J(Q) = S, #S = B = %o, rank 0 input)\n", name, npts, B;
        System(Sprintf("echo 'CLEANUP %o CQ_EXHAUSTED npts=%o' >> ../logs/twcleanup.progress", name, npts));
    else
        printf "CLEANUP curve=%o GAP bound=%o known=%o (extra torsion or untight gcd -> investigate)\n", name, B, #S;
        System(Sprintf("echo 'CLEANUP %o GAP B=%o S=%o' >> ../logs/twcleanup.progress", name, B, #S));
    end if;
end procedure;

// tlabel: [infty-label, "inv"/"mul", c] with T = c/x resp. T = x/c
DoCurve("tw12", F12, [* "T=0", "inv", c12 *]);
DoCurve("tw10", F10, [* "T=infty", "mul", c10 *]);

// (5) unconditional rank bound attempts (LAST; Alarm guards the whole lane)
if DoRank eq 1 then
    printf "RANKPHASE_START (no GRH class-group bounds; Alarm(%o))\n", AlarmS;
    System("echo 'RANKPHASE_START' >> ../logs/twcleanup.progress");
    for entry in [* <"tw12", F12>, <"tw10", F10> *] do
        rb := -1; ok := false;
        try rb := RankBound(Jacobian(HyperellipticCurve(entry[2]))); ok := true; catch e; end try;
        if ok then
            printf "RANKBOUND_UNCOND curve=%o bound=%o\n", entry[1], rb;
            System(Sprintf("echo 'RANKBOUND_UNCOND %o bound=%o' >> ../logs/twcleanup.progress", entry[1], rb));
        else
            printf "RANKBOUND_UNCOND curve=%o FAILED\n", entry[1];
            System(Sprintf("echo 'RANKBOUND_UNCOND %o FAILED' >> ../logs/twcleanup.progress", entry[1]));
        end if;
    end for;
end if;
printf "TWCLEANUP_DONE\n";
System("echo 'TWCLEANUP_DONE' >> ../logs/twcleanup.progress");
quit;
