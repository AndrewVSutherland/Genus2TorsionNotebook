// lane_tw_cleanup4.m — finish the twisted-diagonal cleanup by a class-order
// sieve over finite fields (DivisionPoints is unavailable in genus 3).
// Setting: C = tw12/tw10 odd-degree genus-3 models; J(Q) = J(Q)_tors is a
// finite abelian 2-group (rank 0 UNCONDITIONAL, lane_tw_cleanup.m).  For
// P in C(Q) the class z = [P - infinity] in J(Q) has 2-power order, and for
// every odd good prime p reduction is INJECTIVE on J(Q)_tors, so
// ord(z mod p) = ord(z) exactly.  Sieve: for each p, compute for every
// point Pbar of C(F_p) the order of [Pbar - inf] in J(F_p); keep those of
// pure 2-power order (candidates for reductions of rational points) and
// record the orders.  If across primes the intersection of realized orders
// is contained in {1, 2}, every rational z lies in J(Q)[2]; the reduced
// representative of a factor-subset 2-torsion class has degree = deg of the
// factor product, so the degree<=1 classes are exactly infinity and the
// rational Weierstrass points => C(Q) EXHAUSTED by the known degenerate
// points, unconditionally.
// Usage: cd product/code && magma -b lane_tw_cleanup4.m > ../logs/lane_tw_cleanup4.log
SetColumns(0);
SetMemoryLimit(4*10^9);

QQ := Rationals();
Px<x> := PolynomialRing(QQ);
Pz<T> := PolynomialRing(QQ);

tw12 := -1296*T^8 + 5184*T^7 - 9072*T^6 + 9072*T^5 - 5580*T^4 + 2088*T^3 - 432*T^2 + 36*T;
tw10 := 32*T^7 - 160*T^6 + 256*T^5 - 156*T^4 + 16*T^3 + 16*T^2 - 4*T;

function MonicOdd(f)
    c := LeadingCoefficient(f);
    g := Evaluate(f, Parent(f).1/c) * c^6;
    return Px! [ QQ! Coefficient(g, i) : i in [0..7] ];
end function;

f12flip := Px! [ Coefficient(tw12, 8-i) : i in [0..8] ];
CURVES := [* <"tw12", MonicOdd(f12flip)>, <"tw10", MonicOdd(Px! [ Coefficient(tw10, i) : i in [0..7] ])> *];

NPRIMES := 6;

for ent in CURVES do
    name := ent[1]; f := ent[2];
    printf "\n== %o ==\n", name;
    disc := Integers()! (Numerator(Discriminant(f)) * Denominator(Discriminant(f)));
    C := HyperellipticCurve(f);
    inter := {};    // running intersection of realized pure-2-power orders
    first := true;
    used := [];
    p := 2;
    while #used lt NPRIMES do
        p := NextPrime(p);
        if disc mod p eq 0 then continue; end if;
        okp := true;
        S := {};
        try
            Cp := ChangeRing(C, GF(p));
            Jp := Jacobian(Cp);
            Np := Order(Jp);
            m := Np; v := Valuation(m, 2); modd := m div 2^v;
            pts := Points(Cp);
            for P in pts do
                if P[3] eq 0 then Include(~S, 1); continue; end if;  // infinity: z = 0
                z := Jp! [ P, PointsAtInfinity(Cp)[1] ];
                z2 := modd * z;    // kill odd part
                if not (z2 eq Jp!0) then
                    // pure-2-power candidates only: order of z2 = 2-part
                    o := 1; w := z2;
                    while not (w eq Jp!0) do w := 2*w; o *:= 2; end while;
                    // z has pure 2-power order iff modd*z has full order of z,
                    // i.e. iff odd part of ord(z) is 1 iff z = (something)*z2...
                    // conservative: z is a REDUCTION CANDIDATE only if
                    // ord(z) is a 2-power, i.e. modd*z has the same order as z:
                    // equivalent test: (2^v)*z = 0 exactly when ord | 2^v; then
                    // ord(z) = ord(z2) computed above IF (2-power)*z kills z.
                    if (2^v)*z eq Jp!0 then Include(~S, o); end if;
                else
                    // modd*z = 0: ord(z) | modd (odd): 2-power only if z = 0
                    if z eq Jp!0 then Include(~S, 1); end if;
                end if;
            end for;
        catch e okp := false; end try;
        if not okp then continue; end if;
        Append(~used, p);
        printf "p=%o: #C(F_p)=%o, realized pure-2-power orders %o\n", p, Np, Sort(Setseq(S));
        if first then inter := S; first := false; else inter := inter meet S; end if;
    end while;
    printf "ORDER_SIEVE %o over p in %o: intersection = %o\n", name, used, Sort(Setseq(inter));
    if inter subset {1, 2} then
        fac := [ g[1] : g in Factorization(f) ];
        lin := [ g : g in fac | Degree(g) eq 1 ];
        printf "CQ_FINAL %o: every rational class has order <= 2 => C(Q) = {infinity} + rational Weierstrass points (x = %o) - degenerate points EXHAUST C(Q), UNCONDITIONALLY\n",
            name, [ -Coefficient(l,0)/Coefficient(l,1) : l in lin ];
    else
        printf "CLEANUP4_INCONCLUSIVE %o: orders %o survive the sieve\n", name, Sort(Setseq(inter));
    end if;
end for;
printf "TWCLEANUP4_DONE\n";
quit;
