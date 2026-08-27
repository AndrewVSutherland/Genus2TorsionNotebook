// sibenum.m — Lane B same-member MW descent + lattice enumeration on a T5
// pencil member of M(2,2,2,6).  The member's cover cv (3: y^2=W3, 4: y^2=W4)
// is a genus-1 quartic with rational points at infinity (lc is a square on
// every live member).  Origin = infinity; MinimalModel FIRST (non-minimal
// models caused the lift5.m hang); Generators on the minimal model; then
// enumerate the MW lattice to depth nmax (rank 1) / nm2 box (rank 2) /
// depth 10 box (rank 3), map every point back to u, and exact-test the
// remaining three conditions.  A FULL CANDIDATE (all four squares) is a
// (2,2,2,12) sibling candidate -> jackpot protocol (hitverify.m).
// args (magma -b): rn, rd, cv, nmax, nm2, tmax (CPU-seconds soft cap)
SetClassGroupBounds("GRH");
P<u> := PolynomialRing(Rationals());
q := 4*u^2 - 6*u + 3;
r5 := StringToInteger(rn)/StringToInteger(rd);
CV := StringToInteger(cv);
NMAX := StringToInteger(nmax);
NM2 := StringToInteger(nm2);
TMAX := StringToInteger(tmax);
printf "### sibenum member rho'=%o cover W%o nmax=%o nm2=%o tmax=%o\n",
    r5, CV, NMAX, NM2, TMAX;
W1 := r5*(q*r5-(2*u-1));
W2 := q*r5^2-(4*u^2-4*u+2)*r5+(2*u-1);
W3 := r5*((16*u^4-40*u^3+40*u^2-18*u+3)*r5^3+(-16*u^4+32*u^3-28*u^2+10*u-1)*r5^2
        +(8*u^3-12*u^2+10*u-3)*r5+(-2*u+1));
W4 := (16*u^4-40*u^3+40*u^2-18*u+3)*r5^2+(-16*u^3+28*u^2-18*u+4)*r5+(4*u^2-4*u+1);
G5 := (q*r5-(2*u-1))*(q*r5-1);  // surface-identity class: W1*W2*W3*W4 == G5 mod squares
if CV eq 3 then
    Wc := W3; Wother := [W1, W2, W4];        // [conic, conic, quartic]
elif CV eq 4 then
    Wc := W4; Wother := [W1, W2, W3];
else
    Wc := G5; Wother := [W1, W2, W3, W4];    // cv=5: G-cover, test all four
end if;
known := [-97/48, 133/145, 3637/7105];  // known-curve u representations (verified same G2Invariants)
dc := LCM([Denominator(co) : co in Coefficients(Wc)]);
F := dc^2*Wc;
cont := GCD([Integers()!co : co in Coefficients(F)]);
sq := 1;
for pe in Factorization(cont) do sq *:= pe[1]^(2*(pe[2] div 2)); end for;
F := F div sq;
error if Degree(F) ne 4 or Discriminant(F) eq 0, "degenerate quartic";
C := HyperellipticCurve(F);
ptsI := PointsAtInfinity(C);
if #ptsI gt 0 then
    org := Rep(ptsI);
else
    smp := Points(C : Bound := 100000);
    error if #smp eq 0, "no origin point found";
    org := Rep(smp);
end if;
printf "origin %o on y^2 = %o\n", org, F;
E, mp := EllipticCurve(C, org);
Em, phi := MinimalModel(E);
printf "Em = %o cond %o\n", aInvariants(Em), Conductor(Em);
Tg, mt := TorsionSubgroup(Em);
tor := [mt(t) : t in Tg];
printf "torsion %o\n", Invariants(Tg);
t0 := Cputime();
rlo, rhi := RankBounds(Em);
printf "RankBounds [%o,%o] (%.1o s)\n", rlo, rhi, Cputime(t0);
// NO descent-based Generators (hangs on these curves): harvest small points on
// the quartic cover, push to Em, and lattice-reduce.
t0 := Cputime();
PB := 100000;
smp := Points(C : Bound := PB);
emp := [];
for pt in smp do
    ok := true; Pe := Em!0;
    try
        Pe := phi(mp(pt));
    catch e
        ok := false;
    end try;
    if ok then Append(~emp, Pe); end if;
end for;
printf "cover points H<=%o: %o, pushed to Em: %o (%.1o s)\n", PB, #smp, #emp, Cputime(t0);
t0 := Cputime();
free := ReducedBasis(emp);
printf "ReducedBasis: %o pts, heights %o (%.1o s)\n", #free,
    [RealField(6)!Height(g) : g in free], Cputime(t0);
try
    free := Saturation(free, 11);
    printf "saturated at p<=11: heights %o\n", [RealField(6)!Height(g) : g in free];
catch e
    printf "saturation skipped\n";
end try;
free := [g : g in free | Height(g) gt 0.0001];  // drop torsion injected by Saturation
rk := #free;
if rk lt rlo then
    printf "WARNING: small points give rank %o < lower bound %o -- retry bigger bound\n", rk, rlo;
    smp := Points(C : Bound := 10^7);
    emp := [];
    for pt in smp do
        ok := true; Pe := Em!0;
        try Pe := phi(mp(pt)); catch e ok := false; end try;
        if ok then Append(~emp, Pe); end if;
    end for;
    free := ReducedBasis(emp);
    try free := Saturation(free, 11); catch e ; end try;
    rk := #free;
    printf "after H<=1e7 harvest: rank %o, heights %o\n", rk,
        [RealField(6)!Height(g) : g in free];
end if;
if rk lt rlo then
    printf "WARNING: enumerating a rank-%o sublattice of rank-%o group\n", rk, rlo;
end if;
invphi := Inverse(phi);
invmp := Inverse(mp);

nmap := 0; nfail := 0; ncon := 0; nnear := 0; nfull := 0; nkn := 0;
procedure handle(Pm, ~nmap, ~nfail, ~ncon, ~nnear, ~nfull, ~nkn)
    ok := true; uu := 0;
    try
        Pc := invmp(invphi(Pm));
        if Pc[3] eq 0 then ok := false; else uu := Pc[1]/Pc[3]; end if;
    catch e
        ok := false;
    end try;
    if not ok then nfail +:= 1; return; end if;
    nmap +:= 1;
    vals := [Evaluate(Wo, uu) : Wo in Wother];
    if 0 in vals then return; end if;
    sqs := [IsSquare(v) : v in vals];
    if sqs[1] and sqs[2] then
        ncon +:= 1;
        if &and sqs then
            if uu in known then
                nkn +:= 1;
                printf "KNOWN HIT recovered: u = %o\n", uu;
            else
                nfull +:= 1;
                printf "*** FULL CANDIDATE u = %o *** (member %o cover %o)\n", uu, r5, CV;
            end if;
        else
            nnear +:= 1;
            nd := Ilog(10, Max(Abs(Numerator(uu)), Denominator(uu))) + 1;
            if nd le 40 then
                printf "NEAR3 u = %o (quartic fails)\n", uu;
            else
                printf "NEAR3 u ~ %o digits (quartic fails)\n", nd;
            end if;
        end if;
    end if;
end procedure;

tstart := Cputime();
stopped := false;
if rk eq 0 then
    printf "RANK 0: torsion-only enumeration (%o points)\n", #tor;
    for t in tor do
        handle(t, ~nmap, ~nfail, ~ncon, ~nnear, ~nfull, ~nkn);
    end for;
    if rlo eq 0 and rhi eq 0 and nfull eq 0 then
        printf "CLOSURE: member %o cover W%o has rank 0 => all u with W%o square are torsion images; no sibling hit exists via this cover (GRH class bounds used in descent)\n", r5, CV, CV;
    end if;
elif rk eq 1 then
    PP := free[1];
    Q := Em!0;
    for t in tor do handle(t, ~nmap, ~nfail, ~ncon, ~nnear, ~nfull, ~nkn); end for;
    for n in [1..NMAX] do
        Q +:= PP;
        for t in tor do
            handle(Q+t, ~nmap, ~nfail, ~ncon, ~nnear, ~nfull, ~nkn);
            handle(-Q+t, ~nmap, ~nfail, ~ncon, ~nnear, ~nfull, ~nkn);
        end for;
        if n mod 25 eq 0 then
            printf "  n=%o  mapped=%o fail=%o conpass=%o near3=%o full=%o known=%o  (%.1o s)\n",
                n, nmap, nfail, ncon, nnear, nfull, nkn, Cputime(tstart);
        end if;
        if Cputime(tstart) gt TMAX then
            printf "RESUME: time cap hit at n=%o of %o\n", n, NMAX;
            stopped := true; break;
        end if;
    end for;
elif rk eq 2 then
    // center-out row order so a time cap loses only the far shells
    P1 := free[1]; P2 := free[2];
    Q2s := (-NM2)*P2;
    rows := [0];
    for k in [1..NM2] do Append(~rows, k); Append(~rows, -k); end for;
    for n1 in rows do
        S := n1*P1 + Q2s;
        for n2 in [-NM2..NM2] do
            for t in tor do
                handle(S+t, ~nmap, ~nfail, ~ncon, ~nnear, ~nfull, ~nkn);
            end for;
            S +:= P2;
        end for;
        printf "  row n1=%o  mapped=%o fail=%o conpass=%o near3=%o full=%o known=%o  (%.1o s)\n",
            n1, nmap, nfail, ncon, nnear, nfull, nkn, Cputime(tstart);
        if Cputime(tstart) gt TMAX then
            printf "RESUME: time cap hit after row n1=%o (all |n1| < %o complete)\n", n1, Abs(n1);
            stopped := true; break;
        end if;
    end for;
elif rk eq 3 then
    N3 := Min(NM2, 20);
    P1 := free[1]; P2 := free[2]; P3 := free[3];
    slabs := [0];
    for k in [1..N3] do Append(~slabs, k); Append(~slabs, -k); end for;
    for n1 in slabs do
        R1 := n1*P1;
        for n2 in [-N3..N3] do
            R2 := R1 + n2*P2;
            S := R2 + (-N3)*P3;
            for n3 in [-N3..N3] do
                for t in tor do
                    handle(S+t, ~nmap, ~nfail, ~ncon, ~nnear, ~nfull, ~nkn);
                end for;
                S +:= P3;
            end for;
        end for;
        printf "  slab n1=%o  mapped=%o fail=%o conpass=%o near3=%o full=%o known=%o  (%.1o s)\n",
            n1, nmap, nfail, ncon, nnear, nfull, nkn, Cputime(tstart);
        if Cputime(tstart) gt TMAX then
            printf "RESUME: time cap hit after slab n1=%o (all |n1| < %o complete)\n", n1, Abs(n1);
            stopped := true; break;
        end if;
    end for;
else
    printf "rank %o unhandled\n", rk;
end if;
printf "DONE%o member %o cover W%o: mapped=%o fail=%o conic-pair-pass=%o near3=%o FULL=%o known=%o  total %.1o s\n",
    stopped select " (CAPPED)" else "", r5, CV, nmap, nfail, ncon, nnear, nfull, nkn, Cputime(tstart);
quit;
