//////////////////////////////////////////////////////////////////////
// Stage 2b: point enumeration + candidate factory on the genus-2
// S3-quotient E4 (rank 3), plus split test of its Jacobian.
//
// E4 plane model: mp4Q(Y)=0 over Q(e) (data/contact6_m612_E4_mp4Q.txt),
// hyperelliptic form y^2 = (1/4)x^6-(9/2)x^5+(81/2)x^4-243x^3-729x^2
//                          +8748x+26244  (stage 2), integral twist:
//                y^2 = x^6-18x^5+162x^4-972x^3-2916x^2+34992x+104976.
//
// Pipeline per rational point:
//   (x,y) on H3 -> pull back to plane (Y0,e0)
//   FILTER 1 (C3 lift, free): mp8=mp4(w^2) => lift iff Y0 is a square.
//   FILTER 2 (P8 membership): e0 = e(u) for rational u, via the two
//     nested quadratics t(u)=4(u^2+u-6)/(u^2+6),
//     e(t) = -(25/3)t^2/(t^4-25t^2+1250/3).
// Survivors printed as CANDIDATE e0 for the exact fiber solver.
//
// Also: Jacobian split test via Euler factors at good primes.
//
// Usage: magma -b PtBound:=100000 code/contact6_m612_relative3_E4_points.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
SetMemoryLimit(4*10^9);
if not assigned PtBound then PtBound := 100000;
elif Type(PtBound) eq MonStgElt then PtBound := StringToInteger(PtBound); end if;

Q := Rationals();
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);
mp4Q := Kz ! eval Read("data/contact6_m612_E4_mp4Q.txt");

// plane model over Q
R2<Y,E> := PolynomialRing(Q, 2);
D := LCM([Denominator(Coefficient(mp4Q,i)) : i in [0..4]]);
G := R2!0;
for i in [0..4] do
    c := Coefficient(mp4Q, i) * K!D;
    cn := Numerator(c);
    G +:= (&+[R2 | Coefficient(cn,h)*E^h : h in [0..Degree(cn)]]) * Y^i;
end for;
C := Curve(AffineSpace(R2), G);
bH, H, mpH := IsHyperelliptic(C);
assert bH;
fH, hH := HyperellipticPolynomials(H);
assert hH eq 0;
Px<x> := PolynomialRing(Q);
f4 := Px ! (4*fH);
printf "integral model: y^2 = %o\n", f4;
H3 := HyperellipticCurve(f4);
// H -> H3: (x,y) -> (x,2y)
J3 := Jacobian(H3);
rlo, rhi := RankBounds(J3);
printf "rank bounds on integral model: %o..%o\n", rlo, rhi;
T := TorsionSubgroup(J3);
printf "J torsion: %o\n", Invariants(T);

// split test: Euler factor shape at good primes
print "== split test (chi_p factorization degrees) ==";
disc := Integers()!Discriminant(f4);
for p in [5,7,11,13,17,19,23,29,31,37] do
    if disc mod p eq 0 then continue; end if;
    fp := PolynomialRing(GF(p)) ! f4;
    if not IsSquarefree(fp) then continue; end if;
    Hp := HyperellipticCurve(fp);
    chi := Numerator(ZetaFunction(Hp));
    printf "p=%o : chi factor degrees %o\n", p,
        Sort([Degree(g[1]) : g in Factorization(chi)]);
end for;

// rational points and the candidate pipeline
pts := Points(H3 : Bound := PtBound);
printf "\nH3 rational points (height <= %o): %o\n", PtBound, #pts;

// P8 membership test: exists rational t with e(t)=e0, then rational u with t(u)=t0
function P8Image(e0)
    // e(t) = -(25/3) t^2 / (t^4 - 25 t^2 + 1250/3) = e0
    // <=> e0 t^4 + (-25 e0 + 25/3) t^2 + (1250/3) e0 = 0   (quadratic in t^2)
    if e0 eq 0 then return false, [], "e=0 boundary"; end if;
    Pt<tt> := PolynomialRing(Q);
    quadT2 := e0*tt^2 + (25/3 - 25*e0)*tt + (1250/3)*e0;
    t2roots := [r[1] : r in Roots(quadT2) | r[1] ge 0];
    tvals := [];
    for t2 in t2roots do
        oks, sq := IsSquare(t2);
        if oks then Append(~tvals, sq); Append(~tvals, -sq); end if;
    end for;
    if #tvals eq 0 then return false, [], "no rational t"; end if;
    // t(u) = 4(u^2+u-6)/(u^2+6) = t0  <=>  (4-t0)u^2 + 4u - (24 + 6 t0) = 0
    us := [];
    for t0 in tvals do
        if t0 eq 4 then
            // linear: 4u = 24 + 24 -> u = ... handle
            u0 := (24 + 6*t0)/4; Append(~us, u0); continue;
        end if;
        qu := (4 - t0)*tt^2 + 4*tt - (24 + 6*t0);
        for r in Roots(qu) do Append(~us, r[1]); end for;
    end for;
    if #us eq 0 then return false, [], "no rational u"; end if;
    return true, us, "";
end function;

print "== candidate pipeline ==";
ncand := 0;
for pt in pts do
    // H3 point -> H point -> C point (Y0, e0)
    xx := pt[1]; yy := pt[2]; zzp := pt[3];
    if zzp eq 0 then
        printf "point at infinity %o : (skip; infinite e-fiber boundary)\n", pt;
        continue;
    end if;
    Hpt := H ! [xx, yy/2, 1];
    // pull back through mpH (C -> H birational)
    fib := Hpt @@ mpH;
    rp := RationalPoints(fib);
    if #rp eq 0 then
        printf "point %o : no rational preimage on plane model (bad fiber)\n", pt;
        continue;
    end if;
    for cp in rp do
        Y0 := cp[1]; e0 := cp[2];
        sq, rt := IsSquare(Y0);
        printf "point %o -> (Y,e)=(%o, %o) : Y square? %o", pt, Y0, e0, sq;
        if not sq then printf "  [killed: no C3 lift]\n"; continue; end if;
        okp8, us, why := P8Image(e0);
        if not okp8 then printf "  [killed P8: %o]\n", why; continue; end if;
        ncand +:= 1;
        printf "  *** CANDIDATE e0=%o  u in %o ***\n", e0, us;
    end for;
end for;
printf "\ntotal candidates surviving both filters: %o\n", ncand;
print "E4_POINTS_DONE";
quit;
