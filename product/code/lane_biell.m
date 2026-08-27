// lane_biell.m — bielliptic sweep: for seed elliptic curves E: y^2 = c(x)
// (cubic) and shifts t, the genus-2 curve C_t: y^2 = c_t(x^2) with
// c_t(u) = c(u+t) has J ~ E x E'_t, where E'_t: v^2 = reverse(c_t)(w).
// E'_t roams outside the LMFDB conductor box, so this samples gluing partners
// the seed-pair lane cannot reach.  Funnel only when E'_t has interesting
// torsion.  Usage: magma -b lane_biell.m > ../logs/biell.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/
load "../data/seeds_magma.m";

function CapB(t)
    if t eq [Integers()|2,8] then return <100, 40>;   // <cap, t-height>
    elif t eq [2,6] then return <60, 16>;
    elif t eq [2,4] then return <40, 12>;
    elif t eq [10] then return <100, 24>;
    elif t eq [12] then return <100, 24>;
    elif t eq [9]  then return <100, 20>;
    elif t eq [8]  then return <60, 16>;
    elif t eq [7]  then return <60, 16>;
    else return <0,0>;
    end if;
end function;

function HeightRats(H)
    S := {Rationals()|0};
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Include(~S, a/b); Include(~S, -a/b); end if;
    end for;
    return Sort(Setseq(S));
end function;

counts := AssociativeArray();
nfun := 0; nhit := 0; scanned := 0;
t0 := Cputime();
for s in SEEDS do
    tors := s[4];
    cb := CapB(tors);
    if cb[1] eq 0 then continue; end if;
    key := Sprintf("%o", tors);
    if not IsDefined(counts, key) then counts[key] := 0; end if;
    if counts[key] ge cb[1] then continue; end if;
    counts[key] +:= 1;
    E := EllipticCurve([Rationals()!a : a in s[5]]);
    c := RQx!HyperellipticPolynomials(WeierstrassModel(E));
    oT := #TorsionSubgroup(E);
    minTp := (tors[1] eq 2 and #tors gt 1) select 8 else 3;   // class (a) partners need >= order 8; odd seeds >= 3
    for tv in HeightRats(cb[2]) do
        ct := Evaluate(c, x + tv);
        c0 := Coefficient(ct, 0);
        if c0 eq 0 then continue; end if;
        scanned +:= 1;
        // E'_t: v^2 = reversed cubic  (monicize: x^3+c2 x^2 + a c1 x + a^2 c0
        // for leading coefficient a of the reversed cubic)
        c0R := Coefficient(ct,3); c1R := Coefficient(ct,2);
        c2R := Coefficient(ct,1); aR := Coefficient(ct,0);
        Ep := 0;
        try
            Ep := EllipticCurve([0, c2R, 0, aR*c1R, aR^2*c0R]);
        catch e
            continue;
        end try;
        Tp := TorsionSubgroup(Ep);
        if #Tp lt minTp or oT*#Tp lt 64 then continue; end if;
        oddT := OddPartInvs(Invariants(TorsionSubgroup(E)) cat Invariants(Tp));
        g := Evaluate(ct, x^2);
        den := LCM([Denominator(cc) : cc in Coefficients(g)]);
        g := RQx!(g*den^2);
        if Discriminant(g) eq 0 then continue; end if;
        nfun +:= 1;
        st := Funnel(HyperellipticCurve(g), Sprintf("biell|%o|t=%o|Ep=%o", s[1], tv, Invariants(Tp)) : OddInvs := oddT);
        if st eq "hit" then nhit +:= 1; end if;
    end for;
    printf "PROGRESS biell %o done, scanned %o, funneled %o, hits %o, %o s\n", s[1], scanned, nfun, nhit, Cputime()-t0;
end for;
printf "SEARCH_DONE biell scanned %o funneled %o hits %o %o s\n", scanned, nfun, nhit, Cputime()-t0;
quit;
