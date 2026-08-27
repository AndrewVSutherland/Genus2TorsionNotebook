// lane_q13c3.m — [13]-split via the chi_d-twisted (3,3)-route: scan the
// X1(13) hyperelliptic-fiber quadratic points for E with
//     rho_{E,3} =~ rho_{E,3} (x) chi_d,   d = twist(E^sigma/E),
// the necessary+almost-sufficient rep condition for an anti-3-congruence
// between E and E^sigma = E^{(d)} (the 2-glue is Kani-degenerate for every
// fibral point — campaign note 2026-08-26 — but an anti-isometry on E[3] is
// never twist-induced, so a (3,3)-glue would be nondegenerate).
// Necessary condition scanned here: a_p(E) = 0 mod 3 at EVERY good prime p
// of K with chi_d(p) = -1.  Survivors are printed for the construction
// stage (BHLS 3-glue over K + descent, to be implemented on demand).
// Usage: magma -b SH:=80 Part:=1 NParts:=1 lane_q13c3.m > ../logs/q13c3.log
SetColumns(0);
if not assigned MemGB then MemGB := 8; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned SH then SH := 60; elif Type(SH) eq MonStgElt then SH := StringToInteger(SH); end if;
if not assigned Part then Part := 1; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;

Q := Rationals();
Px<x> := PolynomialRing(Q);
P2<R,S> := PolynomialRing(Q, 2);
F13 := R^3 - R^2*S^4 + 5*R^2*S^3 - 9*R^2*S^2 + 4*R^2*S - 2*R^2 - R*S^3 + 6*R*S^2 - 3*R*S + R - S^3;
A2q<uu,vv> := AffineSpace(Q,2);
hcr := hom< P2 -> CoordinateRing(A2q) | [uu,vv] >;
CX := ProjectiveClosure(Curve(A2q, hcr(F13)));
okh, XH, mph := IsHyperelliptic(CX);
error if not okh, "X1(13) hyperelliptic setup failed";
hpol := Px!HyperellipticPolynomials(XH);
DPm := DefiningPolynomials(mph);
RA := CoordinateRing(Ambient(CX));
hz1 := hom< RA -> P2 | [P2.1, P2.2, 1] >;
XPa := hz1(DPm[1]); ZPa := hz1(DPm[3]);
printf "X1(13) fiber machinery ready\n";

function HeightRats(H)
    S0 := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S0, Q!a/b); Append(~S0, Q!-a/b); end if;
    end for;
    return S0;
end function;

svals := HeightRats(SH);
npts := 0; nsurv := 0;
t0c := Cputime();
for si in [1..#svals] do
    if si mod NParts ne (Part - 1) then continue; end if;
    s0 := svals[si];
    hv := Evaluate(hpol, s0);
    if hv eq 0 or IsSquare(hv) then continue; end if;
    dn := Numerator(hv)*Denominator(hv);
    K := QuadraticField(Squarefree(dn));
    w := K.1;
    R2K := PolynomialRing(K, 2);
    h2 := hom< P2 -> R2K | [R2K.1, R2K.2] >;
    V := [];
    try
        IK := ideal< R2K | h2(F13), h2(XPa) - s0*h2(ZPa) >;
        V := Variety(IK);
    catch e;
        continue;
    end try;
    E := 0; okE := false;
    for pv in V do
        r0 := pv[1]; sr := pv[2];
        if r0 in Q and sr in Q then continue; end if;
        c0 := sr*(r0-1); b0 := r0*c0;
        if b0 eq 0 then continue; end if;
        E0 := 0; okb := true;
        try E0 := EllipticCurve([1-c0, -b0, -b0, 0, 0]); catch e; okb := false; end try;
        if not okb then continue; end if;
        okO := false;
        try okO := Order(E0![0,0]) eq 13; catch e; end try;
        if okO then E := E0; okE := true; break; end if;
    end for;
    if not okE then continue; end if;
    npts +:= 1;
    sig := hom< K -> K | -w >;
    Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
    istw, dtw := IsQuadraticTwist(E, Es);
    if not istw then continue; end if;    // non-fibral anomaly; skip here
    // scan: a_p(E) mod 3 at chi_d-inert primes of K
    OK := Integers(K);
    alive := true; ninert := 0; ngood := 0;
    for p in PrimesInInterval(5, 500) do
        if ninert ge 12 then break; end if;
        dec := Decomposition(OK, p);
        for pr0 in dec do
            pr := pr0[1];
            kp, red := ResidueClassField(pr);
            okr := true; ap := 0; chid := 0;
            try
                dr := red(dtw);
                if dr eq 0 then okr := false; end if;
                if okr then
                    chid := IsSquare(dr) select 1 else -1;
                    ap := TraceOfFrobenius(Reduction(E, pr));
                end if;
            catch e;
                okr := false;
            end try;
            if not okr then continue; end if;
            ngood +:= 1;
            if chid eq -1 then
                ninert +:= 1;
                if ap mod 3 ne 0 then alive := false; break; end if;
            end if;
        end for;
        if not alive then break; end if;
    end for;
    if alive and ninert ge 8 then
        nsurv +:= 1;
        printf "SURV s0=%o D=%o (a_p = 0 mod 3 at %o chi_d-inert primes)\n",
            s0, Squarefree(dn), ninert;
    end if;
    if npts mod 200 eq 0 then
        printf "PROGRESS pts=%o surv=%o %os\n", npts, nsurv, Cputime()-t0c;
    end if;
end for;
printf "SEARCH_DONE q13c3 SH=%o part %o/%o pts=%o survivors=%o %.1o s\n",
    SH, Part, NParts, npts, nsurv, Cputime()-t0c;
quit;
