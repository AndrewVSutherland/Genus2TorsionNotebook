// lane_tw_cleanup2.m — tighten the twisted-diagonal cleanup (2026-08-14):
// lane_tw_cleanup.m established RANK 0 UNCONDITIONALLY for both genus-3
// Jacobians but its 6-prime gcd bounds (128 / 32) exceeded the known
// subgroups (2 / 8).  Here: (1) the EXACT rational 2-torsion from the
// odd-model factorization (odd-degree hyperelliptic: #J(Q)[2] = 2^(k-1),
// k = number of irreducible factors); (2) extended gcd of #J(F_p) over
// good p <= 500 with early exit at the 2-torsion order; (3) if they meet:
// J(Q) = J(Q)[2] exactly (rank 0), every element's reduced Mumford
// representative is a factor-subset divisor, and the degree<=1 elements
// are exactly infinity + the rational Weierstrass points => C(Q) is
// EXHAUSTED by the known degenerate points, unconditionally.
// Usage: cd product/code && magma -b lane_tw_cleanup2.m > ../logs/lane_tw_cleanup2.log
SetColumns(0);
SetMemoryLimit(4*10^9);

QQ := Rationals();
Px<x> := PolynomialRing(QQ);
Pz<T> := PolynomialRing(QQ);

// original twisted kernels (lane_misc2.m:19-20)
tw12 := -1296*T^8 + 5184*T^7 - 9072*T^6 + 9072*T^5 - 5580*T^4 + 2088*T^3 - 432*T^2 + 36*T;
tw10 := 32*T^7 - 160*T^6 + 256*T^5 - 156*T^4 + 16*T^3 + 16*T^2 - 4*T;

// odd-degree monic integral models
// tw12: T | f, flip T -> 1/x: g(x) = x^8 * f(1/x) has degree 7; strip square content
f12flip := Px! [ Coefficient(tw12, 8-i) : i in [0..8] ];   // coefficients reversed
error if Degree(f12flip) ne 7, "tw12 flip not degree 7";
// strip square rational content
cont := &+[ 0 ] + 1;  // placeholder
function MonicOdd(f)
    // scale w^2 = f(T), deg 7, lc c: x = c T, y = c^3 w -> y^2 = c^6 f(x/c) monic
    c := LeadingCoefficient(f);
    g := Evaluate(f, Parent(f).1/c) * c^6;
    g := Px! [ QQ! Coefficient(g, i) : i in [0..7] ];
    error if not IsMonic(g) or Degree(g) ne 7, "monicization failed";
    // clear square denominators: substitute x -> x/s^... (coeffs should be integral already for our cases)
    return g;
end function;

CURVES := [* <"tw12", MonicOdd(f12flip)>, <"tw10", MonicOdd(Px! [ Coefficient(tw10, i) : i in [0..7] ])> *];

for ent in CURVES do
    name := ent[1]; f := ent[2];
    printf "\n== %o: y^2 = %o ==\n", name, f;
    fac := Factorization(f);
    k := #fac;
    printf "FACTORS %o: %o irreducible factors, degrees %o\n", name, k, [ Degree(g[1]) : g in fac ];
    r2 := k - 1;
    printf "EXACT2TORS %o: J(Q)[2] = (Z/2)^%o, order %o\n", name, r2, 2^r2;
    target := 2^r2;
    C := HyperellipticCurve(f);
    disc := Integers()! (Numerator(Discriminant(f)) * Denominator(Discriminant(f)));
    g := 0;
    used := [];
    minv2 := 1000; minv2p := 0;
    for p in [ q : q in [3..500] | IsPrime(q) ] do
        if disc mod p eq 0 then continue; end if;
        okz := true; Np := 0;
        try
            Cp := ChangeRing(C, GF(p));
            Zn := Numerator(ZetaFunction(Cp));
            Np := Integers()! Evaluate(Zn, 1);
        catch e okz := false; end try;
        if not okz or Np eq 0 then continue; end if;
        Append(~used, p);
        g := Gcd(g, Np);
        v2 := Valuation(Np, 2);
        if v2 lt minv2 then minv2 := v2; minv2p := p; end if;
        if g eq target then break; end if;
    end for;
    printf "TORSBOUND2 %o: gcd over %o primes (<= %o) = %o (min v2 = %o at p=%o)\n",
        name, #used, used[#used], g, minv2, minv2p;
    if g eq target then
        printf "TORSION_EXACT %o: J(Q)_tors = J(Q)[2] = (Z/2)^%o (rank 0 UNCONDITIONAL from lane_tw_cleanup.m)\n", name, r2;
        // elements = factor-subset classes; reduced representative of the class
        // of a monic even-subset divisor product g1 (deg <= 3 after reduction
        // g -> f/g when deg > 3); degree-1 representatives <-> rational
        // Weierstrass points
        lin := [ g[1] : g in fac | Degree(g[1]) eq 1 ];
        printf "CQ_EXHAUSTED %o: C(Q) = {infinity} + %o rational Weierstrass points %o - all degenerate\n",
            name, #lin, [ -Coefficient(l,0) : l in lin ];
    else
        printf "TORSION_GAP %o: gcd %o > 2-torsion %o - possible rational 2^k-torsion beyond J[2]; needs halving analysis\n", name, g, target;
    end if;
end for;
printf "TWCLEANUP2_DONE\n";
quit;
