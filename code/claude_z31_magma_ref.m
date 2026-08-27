////////////////////////////////////////////////////////////////////////
// claude_z31_magma_ref.m -- Task B1: Magma reference implementation and
// validation vectors for the Z/31 torsion fiber sieve (charts U' and I).
//
// Chart U' (9 parameters over GF(P)):
//     u = x^2 + u1*x + u0        (monic, deg 2)
//     v = v1*x + v0
//     w = w4*x^4 + w3*x^3 + w2*x^2 + w1*x + w0
//     f = v^2 + u*w              (deg 6 iff w4 != 0; lc(f) = w4)
//   By construction f - v^2 = u*w, so v^2 = f mod u EXACTLY and (u,v)
//   is a valid Mumford pair on y^2 = f: it represents the degree-2
//   affine divisor { (x1, v(x1)), (x2, v(x2)) : u(xi) = 0 } and the
//   Jacobian class DD = [ D_affine - inf+ - inf- ]   (third Mumford
//   component d = 2).
//
// MUMFORD CONSTRUCTION ON EVEN (degree-6) MODELS -- part (a):
//   On even models Magma represents Jacobian points as (a(x), b(x), d).
//   For our class d = 2 = deg(a).  MakeDD below tries, in order (first
//   success wins; per-method counts reported in CHARTU_SUMMARY as
//   makeDD[m1,m2,m3]):
//     m1:  DD := J![u, v];          -- sequence coercion; on even models
//                                      Magma takes d = deg(u) = 2, i.e.
//                                      exactly the class above
//     m2:  DD := elt<J | u, v, 2>;  -- explicit degree-2 constructor
//     m3:  points-difference (P1-I1)+(P2-I2) via roots of u and the two
//          points at infinity (needs u split; diagnostic fallback only)
//   The empirical result on this Magma is recorded in the run logs; see
//   the B1 report.
//
// KUMMER CONVENTION -- part (c):
//   K := KummerSurface(J); PK := K!DD; coords = Eltseq(PK) = [k1..k4].
//   Normalization used in vectors_kummer TSV: first nonzero coord := 1.
//   The script also VERIFIES the Flynn convention
//     (k1:k2:k3:k4) = (1 : x1+x2 : x1*x2 : beta0)
//                   = (1 : -u1 : u0 : (F0 - 2*y1y2)/(u1^2 - 4*u0)),
//   where, with s1 = -u1, s2 = u0 and f = sum_i f_i x^i:
//     F0   = 2*f0 + f1*s1 + 2*f2*s2 + f3*s1*s2 + 2*f4*s2^2
//            + f5*s1*s2^2 + 2*f6*s2^3
//     y1y2 = v1^2*s2 + v1*v0*s1 + v0^2      (= v(x1)*v(x2))
//   Match counters conv_k23 (k2=-u1, k3=u0) and conv_k4 (beta0 formula)
//   are reported in CHARTU_SUMMARY.  K!(J!0) is probed in smoke mode
//   (expected (0:0:0:1)).
//
// Modes (magma -b Mode:=... [PP:=...] claude_z31_magma_ref.m):
//   smoke                quick validation: mini chartU at P=101, Order()
//                        timing probes at P=1009/2003/4003, full part (d),
//                        mini chartI at P=1009
//   chartU PP:=P         parts (a)-(c) for one prime  ->
//                        vectors_chartU_P<P>.tsv:
//                          P u1 u0 v1 v0 w4 w3 w2 w1 w0 order div31flag
//                        vectors_kummer_P<P>.tsv:
//                          same 12 cols + k1 k2 k3 k4 (DD, normalized)
//                          + k1 k2 k3 k4 (31*DD, normalized)
//   rm                   part (d): exact chart-U' tuple of the RM witness
//                        y^2 = F (1830.2.a.q) and its order-31 fibers at
//                        P in {1009,2003,4001,4003}
//   chartI PP:=P         part (e) for one prime -> vectors_chartI_P<P>.tsv:
//                          P c4 c3 c2 c1 c0 N_inf   (N_inf = 0 if > 200)
//                        rows 1-3 are the depressed reductions of the
//                        pell-cf-order validation vectors f14/f18/f28guess
//                        (x -> x - c5/6 kills c5; it fixes both points at
//                        infinity, hence preserves ord(D_inf) exactly);
//                        rows after the 300 random ones are FORCED
//                        small-order rows (0 < N_inf <= 200): CF-filtered
//                        random search, Jacobian-confirmed, so the C
//                        code's positive path is tested on more than the
//                        three knowns.
////////////////////////////////////////////////////////////////////////

if not assigned Mode then Mode := "smoke"; end if;
if not assigned PP then PP := 101; elif Type(PP) eq MonStgElt then PP := StringToInteger(PP); end if;
if not assigned NRand then NRand := 400; elif Type(NRand) eq MonStgElt then NRand := StringToInteger(NRand); end if;
if not assigned NInfRand then NInfRand := 300; elif Type(NInfRand) eq MonStgElt then NInfRand := StringToInteger(NInfRand); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
SetColumns(0);

printf "Z31REF Mode=%o PP=%o NRand=%o NInfRand=%o\n", Mode, PP, NRand, NInfRand;

//////////////////////////////////////////////////////////////// helpers

// Mumford-pair -> Jacobian point on an even (deg-6) model; see header.
function MakeDD(J, u, v)
    DD := J!0; meth := 0; ok := false;
    try
        DD := J![u, v]; ok := true; meth := 1;
    catch e
        ok := false;
    end try;
    if not ok then
        try
            DD := elt<J | u, v, 2>; ok := true; meth := 2;
        catch e
            ok := false;
        end try;
    end if;
    if not ok then
        try
            C := Curve(J);
            rts := Roots(u);
            Pinf := PointsAtInfinity(C);
            if #rts eq 2 and #Pinf eq 2 then
                x1 := rts[1][1]; x2 := rts[2][1];
                P1 := C![x1, Evaluate(v, x1)];
                P2 := C![x2, Evaluate(v, x2)];
                DD := (P1 - Pinf[1]) + (P2 - Pinf[2]);
                ok := true; meth := 3;
            end if;
        catch e
            ok := false;
        end try;
    end if;
    if ok then return true, DD, meth; end if;
    return false, _, 0;
end function;

// projective normalization: first nonzero coordinate -> 1
function NormProj(s)
    i := 1;
    while i le #s and s[i] eq 0 do i +:= 1; end while;
    if i gt #s then return s; end if;
    c := s[i];
    return [ t/c : t in s ];
end function;

function KumCoords(pt)
    try
        s := Eltseq(pt);
        return s;
    catch e
        return [ pt[i] : i in [1..4] ];
    end try;
end function;

// polynomial part of sqrt(f) for monic sextic f (any field, char != 2)
function SqrtPolyPartFp(f)
    Rp := Parent(f); x := Rp.1;
    s := x^3;
    for k in [1..3] do
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        s := s + (Coefficient(d, 6-k)/(2*Coefficient(s, 3)))*x^(3-k);
    end for;
    return s;
end function;

// exact ord(D_inf) via polynomial CF of sqrt(f), early-abort above cap.
// Port of the validated code/agent_a2_24_cf.m (same recursion, same
// quasi-period stop, INCLUDING deg(a_0)=3 in the total), specialized to
// finite prime fields where every nonzero lc is a unit.
function CFOrderFp(f, cap)
    s := SqrtPolyPartFp(f);
    Rp := Parent(f);
    Pi := Rp!0; Qi := Rp!1; total := 0;
    for i in [0..cap+10] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + s) div Qi;
        dai := Degree(ai);
        if dai lt 0 then return 0; end if;
        if i ge 1 and dai lt 1 then return 0; end if;   // degenerate orbit
        total +:= dai;
        if total gt cap then return 0; end if;          // early abort
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;  // nonexact
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;

///////////////////////////////////////////////////////// chart U' bulk

procedure DoChartU(PP, NRand)
    Z := Integers();
    Fp := GF(PP);
    R<x> := PolynomialRing(Fp);
    SetSeed(31000 + PP);
    fU := Sprintf("vectors_chartU_P%o.tsv", PP);
    fK := Sprintf("vectors_kummer_P%o.tsv", PP);
    System("rm -f " cat fU);
    System("rm -f " cat fK);
    nraw := 0; nvalid := 0; ndiv := 0; nfail := 0; nskipdeg := 0;
    m1 := 0; m2 := 0; m3 := 0;
    krows := 0; kdiv := 0; kfail := 0;
    conv23ok := 0; conv4ok := 0; convtot := 0;
    it := 0; maxit := 200*NRand + 30000;
    t0 := Cputime();
    while true do
        it +:= 1;
        if it gt maxit then printf "CHARTU_WARN P=%o hit maxit=%o with ndiv=%o\n", PP, maxit, ndiv; break; end if;
        inRandomPhase := nraw lt NRand;
        if (not inRandomPhase) and ndiv ge 5 then break; end if;
        u1 := Random(Fp); u0 := Random(Fp); v1 := Random(Fp); v0 := Random(Fp);
        w4 := Random(Fp); w3 := Random(Fp); w2 := Random(Fp); w1 := Random(Fp); w0 := Random(Fp);
        if inRandomPhase then nraw +:= 1; end if;
        u := x^2 + u1*x + u0;
        v := v1*x + v0;
        w := w4*x^4 + w3*x^3 + w2*x^2 + w1*x + w0;
        f := v^2 + u*w;
        if Degree(f) ne 6 or Degree(Gcd(f, Derivative(f))) ne 0 then
            nskipdeg +:= 1; continue;
        end if;
        okc := true;
        try
            C := HyperellipticCurve(f); J := Jacobian(C);
        catch e
            okc := false;
        end try;
        if not okc then nskipdeg +:= 1; continue; end if;
        okD, DD, meth := MakeDD(J, u, v);
        if not okD then nfail +:= 1; continue; end if;
        if meth eq 1 then m1 +:= 1; elif meth eq 2 then m2 +:= 1; else m3 +:= 1; end if;
        N := Order(DD);
        flag := (N mod 31 eq 0) select 1 else 0;
        if (not inRandomPhase) and flag eq 0 then continue; end if;  // top-up records hits only
        nvalid +:= 1;
        if flag eq 1 then ndiv +:= 1; end if;
        row := Sprintf("%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o",
            PP, Z!u1, Z!u0, Z!v1, Z!v0, Z!w4, Z!w3, Z!w2, Z!w1, Z!w0, N, flag);
        fprintf fU, "%o\n", row;
        wantK := (nvalid le 10) or (flag eq 1 and kdiv lt 3);
        if wantK then
            okk := true;
            kD := [Fp|0,0,0,0]; k3 := [Fp|0,0,0,0];
            try
                K := KummerSurface(J);
                kD := NormProj(KumCoords(K!DD));
                S := 31*DD;
                if S eq J!0 then
                    k3 := [Fp|0, 0, 0, 1];
                else
                    k3 := NormProj(KumCoords(K!S));
                end if;
            catch e
                okk := false;
            end try;
            if not okk then
                kfail +:= 1;
            else
                krows +:= 1;
                if flag eq 1 then kdiv +:= 1; end if;
                // convention checks (Flynn; see header)
                convtot +:= 1;
                if kD[1] eq 1 and kD[2] eq -u1 and kD[3] eq u0 then conv23ok +:= 1; end if;
                s1 := -u1; s2 := u0;
                fc := [Coefficient(f, i) : i in [0..6]];
                F0 := 2*fc[1] + fc[2]*s1 + 2*fc[3]*s2 + fc[4]*s1*s2
                      + 2*fc[5]*s2^2 + fc[6]*s1*s2^2 + 2*fc[7]*s2^3;
                y1y2 := v1^2*s2 + v1*v0*s1 + v0^2;
                den := s1^2 - 4*s2;
                if den ne 0 and kD[1] eq 1 then
                    k4pred := (F0 - 2*y1y2)/den;
                    if kD[4] eq k4pred then conv4ok +:= 1; end if;
                end if;
                if krows le 2 then
                    printf "KUMMER_SAMPLE P=%o u1=%o u0=%o normK(DD)=%o normK(31DD)=%o\n",
                        PP, Z!u1, Z!u0, kD, k3;
                end if;
                krow := row cat Sprintf("\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o",
                    Z!kD[1], Z!kD[2], Z!kD[3], Z!kD[4],
                    Z!k3[1], Z!k3[2], Z!k3[3], Z!k3[4]);
                fprintf fK, "%o\n", krow;
            end if;
        end if;
        if nvalid mod 50 eq 0 then
            printf "PROGRESS chartU P=%o tried=%o valid=%o div31=%o t=%os\n", PP, it, nvalid, ndiv, Cputime(t0);
        end if;
    end while;
    printf "CHARTU_SUMMARY P=%o tried=%o raw=%o valid_rows=%o div31=%o skipdeg=%o makeDD[m1,m2,m3]=[%o,%o,%o] makeDD_fails=%o kummer_rows=%o kummer_fails=%o conv_k23=%o/%o conv_k4=%o/%o time=%os\n",
        PP, it, nraw, nvalid, ndiv, nskipdeg, m1, m2, m3, nfail,
        krows, kfail, conv23ok, convtot, conv4ok, convtot, Cputime(t0);
end procedure;

/////////////////////////////////////////////////////////// part (d) RM

procedure DoRM(~dummy)
    Z := Integers();
    QQ := Rationals();
    RQ<X> := PolynomialRing(QQ);
    F := -3356*X^6 + 11364*X^5 - 18347*X^4 + 17202*X^3 - 9863*X^2 + 3264*X - 504;
    uQ := X^2 - 21/20*X + 9/20;
    vQ := 101/200*X + 51/200;
    num := F - vQ^2;
    assert num mod uQ eq 0;
    wQ := num div uQ;
    assert uQ*wQ + vQ^2 eq F;
    tupQ := [Coefficient(uQ, 1), Coefficient(uQ, 0),
             Coefficient(vQ, 1), Coefficient(vQ, 0)]
            cat [Coefficient(wQ, i) : i in [4, 3, 2, 1, 0]];
    printf "RM_CHARTU_TUPLE (u1,u0,v1,v0,w4,w3,w2,w1,w0) = %o\n", tupQ;
    printf "RM_W_POLY w_RM = %o\n", wQ;
    DiscF := Discriminant(F);
    printf "RM_DISC disc(F) = %o\n", DiscF;
    denlcm := LCM([Denominator(c) : c in tupQ]);
    printf "RM_DEN_LCM = %o\n", denlcm;
    for P in [1009, 2003, 4001, 4003] do
        if denlcm mod P eq 0 or Numerator(DiscF) mod P eq 0
           or (Z!Coefficient(F, 6)) mod P eq 0 then
            printf "RMFIBER P=%o SKIPPED (divides denominators/disc/lc)\n", P;
            continue;
        end if;
        Fp := GF(P);
        Rp<xp> := PolynomialRing(Fp);
        fP := Rp![Fp!Coefficient(F, i) : i in [0..6]];
        uP := Rp![Fp!Coefficient(uQ, i) : i in [0..2]];
        vP := Rp![Fp!Coefficient(vQ, i) : i in [0..1]];
        wP := Rp![Fp!Coefficient(wQ, i) : i in [0..4]];
        assert vP^2 + uP*wP eq fP;
        assert Degree(Gcd(fP, Derivative(fP))) eq 0;
        J := Jacobian(HyperellipticCurve(fP));
        ok, DD, meth := MakeDD(J, uP, vP);
        assert ok;
        N := Order(DD);
        printf "RMFIBER P=%o tuple=%o order=%o meth=%o %o\n", P,
            [Z!(Fp!c) : c in tupQ], N, meth,
            (N eq 31) select "PASS" else "FAIL";
    end for;
    printf "RM_DONE\n";
end procedure;

//////////////////////////////////////////////////////// part (e) chart I

procedure DoChartI(PP, NInfRand)
    Z := Integers();
    Fp := GF(PP);
    R<x> := PolynomialRing(Fp);
    SetSeed(77000 + PP);
    fI := Sprintf("vectors_chartI_P%o.tsv", PP);
    System("rm -f " cat fI);
    t0 := Cputime();
    // (1) the three validated pell-cf-order Q-vectors, reduced mod PP and
    //     depressed x -> x - c5/6 (fixes inf+/-, preserves ord(D_inf)).
    knowns := [
        <"f14", (x^2 + 1)*(x^4 + 5*x^2 + 4*x + 4), 14>,
        <"f18", (x^2 - x + 1)*(x^4 - x^3 + 9*x^2 + 8*x - 8), 18>,
        <"f28guess", x^6 + 2*x^5 - 5*x^4 - 14*x^3 - 3*x^2 + 24*x + 28, 7>
    ];
    npass := 0;
    for kv in knowns do
        g := kv[2];
        c5 := Coefficient(g, 5);
        g2 := Evaluate(g, x - c5/6);
        assert Coefficient(g2, 5) eq 0;
        assert Degree(Gcd(g2, Derivative(g2))) eq 0;
        C := HyperellipticCurve(g2);
        Pinf := PointsAtInfinity(C);
        assert #Pinf eq 2;
        DD := Pinf[1] - Pinf[2];
        N := Order(DD);
        cf := CFOrderFp(g2, 200);
        okk := (N eq kv[3]) and (cf eq kv[3]);
        if okk then npass +:= 1; end if;
        printf "CHARTI_KNOWN P=%o %o depressed_c=[%o,%o,%o,%o,%o] jac_order=%o cf_order=%o expected=%o %o\n",
            PP, kv[1], Z!Coefficient(g2, 4), Z!Coefficient(g2, 3), Z!Coefficient(g2, 2),
            Z!Coefficient(g2, 1), Z!Coefficient(g2, 0), N, cf, kv[3],
            okk select "PASS" else "FAIL";
        Ncap := (N gt 200) select 0 else N;
        fprintf fI, "%o\t%o\t%o\t%o\t%o\t%o\t%o\n", PP,
            Z!Coefficient(g2, 4), Z!Coefficient(g2, 3), Z!Coefficient(g2, 2),
            Z!Coefficient(g2, 1), Z!Coefficient(g2, 0), Ncap;
    end for;
    // (2) random monic depressed sextics
    nrows := 0; nskip := 0; agree := 0; mism := 0; nsmall := 0;
    for iter in [1..NInfRand] do
        c4 := Random(Fp); c3 := Random(Fp); c2 := Random(Fp);
        c1 := Random(Fp); c0 := Random(Fp);
        f := x^6 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0;
        if Degree(Gcd(f, Derivative(f))) ne 0 then nskip +:= 1; continue; end if;
        okc := true;
        try
            C := HyperellipticCurve(f);
            Pinf := PointsAtInfinity(C);
            J := Jacobian(C);
            DD := Pinf[1] - Pinf[2];
        catch e
            okc := false;
        end try;
        if not okc or #Pinf ne 2 then nskip +:= 1; continue; end if;
        N := Order(DD);
        Ncap := (N gt 200) select 0 else N;
        if Ncap ne 0 then nsmall +:= 1; end if;
        cf := CFOrderFp(f, 200);
        if cf eq Ncap then
            agree +:= 1;
        else
            mism +:= 1;
            printf "CF_MISMATCH P=%o c=[%o,%o,%o,%o,%o] jac=%o cap=%o cf=%o\n",
                PP, Z!c4, Z!c3, Z!c2, Z!c1, Z!c0, N, Ncap, cf;
        end if;
        fprintf fI, "%o\t%o\t%o\t%o\t%o\t%o\t%o\n", PP,
            Z!c4, Z!c3, Z!c2, Z!c1, Z!c0, Ncap;
        nrows +:= 1;
        if nrows mod 50 eq 0 then
            printf "PROGRESS chartI P=%o rows=%o t=%os\n", PP, nrows, Cputime(t0);
        end if;
    end for;
    // (3) FORCED small-order rows: search until >= 10 rows with
    //     0 < N_inf <= 200.  CF (cheap, exact) is the filter; every hit is
    //     confirmed by Jacobian arithmetic before being written.
    nforce := 0; ftries := 0; fdisagree := 0;
    while nforce lt 10 and ftries lt 400000 do
        ftries +:= 1;
        if ftries mod 25000 eq 0 then
            printf "PROGRESS chartI-force P=%o tries=%o forced=%o t=%os\n", PP, ftries, nforce, Cputime(t0);
        end if;
        c4 := Random(Fp); c3 := Random(Fp); c2 := Random(Fp);
        c1 := Random(Fp); c0 := Random(Fp);
        f := x^6 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0;
        if Degree(Gcd(f, Derivative(f))) ne 0 then continue; end if;
        cf := CFOrderFp(f, 200);
        if cf eq 0 then continue; end if;
        okc := true;
        try
            C := HyperellipticCurve(f);
            Pinf := PointsAtInfinity(C);
            J := Jacobian(C);
            DD := Pinf[1] - Pinf[2];
        catch e
            okc := false;
        end try;
        if not okc or #Pinf ne 2 then continue; end if;
        N := Order(DD);
        if N ne cf then
            fdisagree +:= 1;
            printf "CF_FORCE_MISMATCH P=%o c=[%o,%o,%o,%o,%o] cf=%o jac=%o\n",
                PP, Z!c4, Z!c3, Z!c2, Z!c1, Z!c0, cf, N;
            continue;
        end if;
        fprintf fI, "%o\t%o\t%o\t%o\t%o\t%o\t%o\n", PP,
            Z!c4, Z!c3, Z!c2, Z!c1, Z!c0, N;
        nforce +:= 1;
        printf "CHARTI_FORCED P=%o c=[%o,%o,%o,%o,%o] N_inf=%o tries=%o t=%os\n",
            PP, Z!c4, Z!c3, Z!c2, Z!c1, Z!c0, N, ftries, Cputime(t0);
    end while;
    printf "CHARTI_SUMMARY P=%o knowns_pass=%o/3 raw=%o rows=%o skip=%o small_order=%o cf_agree=%o cf_mismatch=%o forced=%o forced_tries=%o forced_disagree=%o time=%os\n",
        PP, npass, NInfRand, nrows, nskip, nsmall, agree, mism, nforce, ftries, fdisagree, Cputime(t0);
end procedure;

//////////////////////////////////////////////////////// timing probe

procedure TimeOrders(PP, n)
    Fp := GF(PP);
    R<x> := PolynomialRing(Fp);
    SetSeed(999);
    t0 := Cputime(); cnt := 0; tries := 0;
    while cnt lt n and tries lt 100*n do
        tries +:= 1;
        u1 := Random(Fp); u0 := Random(Fp); v1 := Random(Fp); v0 := Random(Fp);
        w4 := Random(Fp); w3 := Random(Fp); w2 := Random(Fp); w1 := Random(Fp); w0 := Random(Fp);
        u := x^2 + u1*x + u0; v := v1*x + v0;
        w := w4*x^4 + w3*x^3 + w2*x^2 + w1*x + w0;
        f := v^2 + u*w;
        if Degree(f) ne 6 or Degree(Gcd(f, Derivative(f))) ne 0 then continue; end if;
        okc := true;
        try
            J := Jacobian(HyperellipticCurve(f));
        catch e
            okc := false;
        end try;
        if not okc then continue; end if;
        okD, DD, meth := MakeDD(J, u, v);
        if not okD then continue; end if;
        N := Order(DD);
        cnt +:= 1;
    end while;
    printf "TIMING P=%o n=%o total=%os per_order=%os\n", PP, cnt, Cputime(t0), Cputime(t0)/Max(cnt, 1);
end procedure;

//////////////////////////////////////////////////////////// smoke bits

procedure SmokeKummerZero(~dummy)
    Fp := GF(101);
    R<x> := PolynomialRing(Fp);
    // fixed nondegenerate tuple
    u := x^2 + 3*x + 5; v := 2*x + 7;
    w := 11*x^4 + 13*x^3 + 17*x^2 + 19*x + 23;
    f := v^2 + u*w;
    assert Degree(f) eq 6 and Degree(Gcd(f, Derivative(f))) eq 0;
    J := Jacobian(HyperellipticCurve(f));
    K := KummerSurface(J);
    try
        printf "K_ZERO_TEST K!(J!0) coords = %o\n", KumCoords(K!(J!0));
    catch e
        printf "K_ZERO_TEST threw\n";
    end try;
    printf "K_EQUATION P=101 sample: %o\n", DefiningPolynomial(K);
end procedure;

/////////////////////////////////////////////////////////////// dispatch

dummy := 0;
if Mode eq "smoke" then
    DoChartU(101, 20);
    SmokeKummerZero(~dummy);
    TimeOrders(1009, 10);
    TimeOrders(2003, 10);
    TimeOrders(4003, 5);
    DoRM(~dummy);
    DoChartI(1009, 10);
elif Mode eq "chartU" then
    DoChartU(PP, NRand);
elif Mode eq "rm" then
    DoRM(~dummy);
elif Mode eq "chartI" then
    DoChartI(PP, NInfRand);
else
    printf "UNKNOWN Mode %o\n", Mode;
end if;
print "ALLDONE";
quit;
