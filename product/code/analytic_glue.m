// analytic_glue.m — analytic (N,N)-gluing of two elliptic curves along an
// anti-isometry of their N-torsion, WITHOUT explicit modular formulas:
//   Lambda' = Lambda_{ExF} + graph_M/N  for M in GL2(Z/N);
//   form N*E_prod on Lambda' (integral+unimodular iff the graph is
//   Lagrangian, i.e. M is an anti-isometry);
//   Frobenius symplectic basis -> small period matrix tau -> Rosenhain ->
//   Igusa-Clebsch -> rationality detection (only Galois-equivariant M's give
//   rational invariants) -> curve over Q up to twist -> twist fixed by
//   matching #C(F_p) = p+1-a_E(p)-a_F(p).
// The output curves C satisfy Jac(C) = (ExF)/graph and hence
// C's Jacobian has rational torsion containing T1 x T2 whenever
// gcd(N, #T1*#T2) = 1.  Everything is certified afterwards by exact
// point-count consistency and TorsionSubgroup (the analytic step needs no
// proof).
// Load from product/code/:  load "analytic_glue.m";
// Main entry:  GlueCandidates(E, F, N : prec := 120, TwistPrimes := 60)

AGprec := 120;

function EProdPair(x, y)
    // product symplectic form in coordinates w.r.t. (w1,w2,w1',w2')
    return x[1]*y[2]-x[2]*y[1] + x[3]*y[4]-x[4]*y[3];
end function;

function NormPeriods(E, CC)
    per := Periods(E : Precision := Precision(CC));
    w1 := CC!per[1]; w2 := CC!per[2];
    if Im(w2/w1) lt 0 then tmp := w1; w1 := w2; w2 := tmp; end if;
    return [w1, w2];
end function;

function RatApprox(z, eps, hbound)
    // z complex; returns ok, rational
    // WARNING (2026-08-13 lesson, see notes/claude_split_2266_1212_session_2026_08_13.md
    // sec. 4.7): with eps ~ 1/hbound this WEAK test accepts essentially every
    // real number (BestApproximation beats 1/hbound generically).  Kept only
    // for reference; all recognition below goes through StrictRatApprox.
    if Abs(Im(z)) gt eps then return false, 0; end if;
    r := Re(z);
    q := BestApproximation(r, hbound);
    if Abs(r - q) gt eps then return false, 0; end if;
    return true, q;
end function;

function StrictRatApprox(z, prec)
    // strict height-relative rational recognition: accept q only when the
    // error beats 10^-(2*height(q)+15), i.e. the match is twice as good as
    // the information content of q — which a random real never achieves.
    // Requires prec >= 2*height + 35; the height ladder stops accordingly.
    // (Same convention as glue_window.m / lane_sigma_strictglue.m.)
    if Abs(Im(z)) gt RealField(prec)!10.0^(-(prec div 3)) then return false, 0; end if;
    r := Re(z);
    for hb in [20, 40, 60, 90, 130, 180, 250, 350, 500] do
        if 2*hb + 20 gt prec then break; end if;
        q := BestApproximation(r, 10^hb);
        hq := Max(Ilog(10, 1+Abs(Numerator(q))), Ilog(10, 1+Abs(Denominator(q))));
        if Abs(r - q) lt RealField(prec)!10.0^(-(2*hq + 15)) then return true, q; end if;
    end for;
    return false, 0;
end function;

function ReduceIC(IC)
    // weighted content reduction of a rational Igusa-Clebsch tuple
    // (weights 2,4,6,10) by primes found via trial division; avoids any
    // full factorization of large invariants.
    wts := [2,4,6,10];
    supp := {};
    for q in IC do
        if q eq 0 then continue; end if;
        for pe in TrialDivision(Numerator(q), 10^6) do Include(~supp, pe[1]); end for;
        for pe in TrialDivision(Denominator(q), 10^6) do Include(~supp, pe[1]); end for;
    end for;
    out := IC;
    for p in supp do
        e := Minimum([ Floor(Valuation(out[i], p)/wts[i]) : i in [1..4] | out[i] ne 0 ]);
        if e ne 0 then
            out := [ out[i] eq 0 select 0 else out[i]/p^(e*wts[i]) : i in [1..4] ];
        end if;
    end for;
    return out;
end function;

function GlueScan(E, F, N : prec := 60, Verbose := true)
    // FAST stage: enumerate Lagrangian M, detect Galois-equivariant gluings
    // by rational Igusa-Clebsch invariants at moderate precision; NO curve
    // construction.  Returns list of < [m11,m12,m21,m22], q1, q2, q3 >.
    CC := ComplexField(prec);
    wE := NormPeriods(E, CC); wF := NormPeriods(F, CC);
    ZZ := Integers(); QQ := Rationals();
    out := [* *];
    for m11, m12, m21, m22 in [0..N-1] do
        if (m11*m22 - m12*m21 + 1) mod N ne 0 then continue; end if;
        M := Matrix(Integers(N), 2, 2, [m11,m12,m21,m22]);
        if not IsInvertible(M) then continue; end if;
        rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
                  [QQ|1/N, 0, m11/N, m12/N], [QQ|0, 1/N, m21/N, m22/N] ];
        H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| N*c : c in r ] : r in rows ]));
        B := [ [ QQ | H[i][j]/N : j in [1..4] ] : i in [1..4] ];
        J := Matrix(QQ, 4,4, [ [ N*EProdPair(B[i], B[j]) : j in [1..4] ] : i in [1..4] ]);
        if not forall{ <i,j> : i,j in [1..4] | Denominator(J[i][j]) eq 1 } then continue; end if;
        JZ := Matrix(ZZ, 4,4, [ [ ZZ!J[i][j] : j in [1..4] ] : i in [1..4] ]);
        if Abs(Determinant(JZ)) ne 1 then continue; end if;
        F0, T := FrobeniusFormAlternating(JZ);
        if F0[1][3] ne 1 or F0[2][4] ne 1 then continue; end if;
        Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
        function colv(v)
            return [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ];
        end function;
        cols := [ colv(Cb[i]) : i in [1..4] ];
        PA := Matrix(CC, 2,2, [ cols[1][1], cols[2][1], cols[1][2], cols[2][2] ]);
        PB := Matrix(CC, 2,2, [ cols[3][1], cols[4][1], cols[3][2], cols[4][2] ]);
        okt := true; tau := 0;
        try tau := PA^-1*PB; catch e okt := false; end try;
        if not okt then continue; end if;
        function imposdef(t)
            a := Im(t[1][1]); d := Im(t[2][2]); b := Im(t[1][2]);
            return a gt 0 and a*d - b^2 gt 0;
        end function;
        if not imposdef(tau) then
            try tau := PB^-1*PA; catch e okt := false; end try;
            if not okt or not imposdef(tau) then continue; end if;
        end if;
        okr := true; ros := [];
        try ros := RosenhainInvariants(tau); catch e okr := false; end try;
        if not okr then continue; end if;
        if exists{ r : r in ros | Abs(r) lt 10.0^(-8) or Abs(r-1) lt 10.0^(-8) } then continue; end if;
        PCx<xx> := PolynomialRing(CC);
        g := xx*(xx-1)*&*[ xx - r : r in ros ];
        IC := IgusaClebschInvariants(HyperellipticCurve(g));
        if Abs(IC[4]) lt 10.0^(-prec div 2) then continue; end if;
        j1 := IC[1]^5/IC[4]; j2 := IC[1]^3*IC[2]/IC[4]; j3 := IC[1]^2*IC[3]/IC[4];
        ok1, q1 := StrictRatApprox(j1, prec);
        if not ok1 then continue; end if;
        ok2, q2 := StrictRatApprox(j2, prec);
        ok3, q3 := StrictRatApprox(j3, prec);
        if not (ok2 and ok3) then continue; end if;
        if Verbose then printf "GLUESCAN M=[%o,%o;%o,%o] heights %o %o %o\n", m11,m12,m21,m22,
            Ilog(10, 1+Abs(Numerator(q1))), Ilog(10, 1+Abs(Numerator(q2))), Ilog(10, 1+Abs(Numerator(q3))); end if;
        Append(~out, < [m11,m12,m21,m22], q1, q2, q3 >);
    end for;
    return out;
end function;

function GlueCandidates(E, F, N : prec := 120, TwistPrimes := 60, Verbose := true)
    // returns list of <curve over Q, M> ; curves deduped by G2 invariants
    CC := ComplexField(prec);
    RQxL := PolynomialRing(Rationals());
    wE := NormPeriods(E, CC); wF := NormPeriods(F, CC);
    ZZ := Integers(); QQ := Rationals();
    out := [* *];
    seen := {};
    aE := AssociativeArray(); aF := AssociativeArray();
    NE := Conductor(E); NF := Conductor(F);
    for p in PrimesUpTo(TwistPrimes) do
        if (NE*NF) mod p ne 0 then
            aE[p] := TraceOfFrobenius(E, p); aF[p] := TraceOfFrobenius(F, p);
        end if;
    end for;
    // squarefree twist candidates from 2*rad(NE*NF)
    rad := &*[ZZ| q : q in PrimeDivisors(2*NE*NF) ];
    tws := [ ZZ | d : d in Divisors(rad) ];
    tws := tws cat [ -d : d in tws ];
    nM := 0; nlag := 0; nrat := 0;
    for m11, m12, m21, m22 in [0..N-1] do
        // Lagrangian graph <=> det = -1 mod N (then N*E_prod is integral
        // unimodular on Lambda'); cheap precheck before any lattice work
        if (m11*m22 - m12*m21 + 1) mod N ne 0 then continue; end if;
        M := Matrix(Integers(N), 2, 2, [m11,m12,m21,m22]);
        if not IsInvertible(M) then continue; end if;
        nM +:= 1;
        // lattice Lambda' in rational coords: rows
        rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
                  [QQ|1/N, 0, m11/N, m12/N], [QQ|0, 1/N, m21/N, m22/N] ];
        H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| N*c : c in r ] : r in rows ]));
        B := [ [ QQ | H[i][j]/N : j in [1..4] ] : i in [1..4] ];
        // form N*E_prod on the basis
        J := Matrix(QQ, 4,4, [ [ N*EProdPair(B[i], B[j]) : j in [1..4] ] : i in [1..4] ]);
        if not forall{ <i,j> : i,j in [1..4] | Denominator(J[i][j]) eq 1 } then continue; end if;
        JZ := Matrix(ZZ, 4,4, [ [ ZZ!J[i][j] : j in [1..4] ] : i in [1..4] ]);
        if Abs(Determinant(JZ)) ne 1 then continue; end if;
        nlag +:= 1;
        F0, T := FrobeniusFormAlternating(JZ);
        error if F0[1][3] ne 1 or F0[2][4] ne 1, "nonprincipal Frobenius form";
        // new basis c_i = sum_j T[i][j] B[j]; F0 should be [[0,I],[-I,0]]
        Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
        // period matrix columns
        function colv(v)
            return [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ];
        end function;
        cols := [ colv(Cb[i]) : i in [1..4] ];
        PA := Matrix(CC, 2,2, [ cols[1][1], cols[2][1], cols[1][2], cols[2][2] ]);
        PB := Matrix(CC, 2,2, [ cols[3][1], cols[4][1], cols[3][2], cols[4][2] ]);
        okt := true; tau := 0;
        try tau := PA^-1*PB; catch e okt := false; end try;
        if not okt then continue; end if;
        // symmetrize check and positivity; fall back to swapped order
        function imposdef(t)
            a := Im(t[1][1]); d := Im(t[2][2]); b := Im(t[1][2]);
            return a gt 0 and a*d - b^2 gt 0;
        end function;
        if not imposdef(tau) then
            try tau := PB^-1*PA; catch e okt := false; end try;
            if not okt or not imposdef(tau) then continue; end if;
        end if;
        okr := true; ros := [];
        try ros := RosenhainInvariants(tau); catch e okr := false; end try;
        if not okr then continue; end if;
        PCx<xx> := PolynomialRing(CC);
        g := xx*(xx-1)*&*[ xx - r : r in ros ];
        if exists{ r : r in ros | Abs(r) lt 10.0^(-10) or Abs(r-1) lt 10.0^(-10) } then continue; end if;
        IC := IgusaClebschInvariants(HyperellipticCurve(g));
        if Abs(IC[4]) lt 10.0^(-prec div 2) then continue; end if;  // degenerate/product
        j1 := IC[1]^5/IC[4]; j2 := IC[1]^3*IC[2]/IC[4]; j3 := IC[1]^2*IC[3]/IC[4];
        ok1, q1 := StrictRatApprox(j1, prec);
        if not ok1 then continue; end if;
        ok2, q2 := StrictRatApprox(j2, prec);
        ok3, q3 := StrictRatApprox(j3, prec);
        if not (ok2 and ok3) then continue; end if;
        nrat +:= 1;
        if Verbose then printf "GLUE M=[%o,%o;%o,%o] rational invariants %o %o %o\n", m11,m12,m21,m22, q1, q2, q3; end if;
        if q1 eq 0 then printf "GLUE: I2 = 0 case, skipping (needs other normalization)\n"; continue; end if;
        // normalize I2 = 1: tuple [1, I4/I2^2, I6/I2^3, I10/I2^5] carries the
        // honest invariant heights (the earlier [q1,...,q1^4] form inflated
        // heights ~4x via unremovable large-prime content)
        ICq := ReduceIC([ QQ | 1, q2/q1, q3/q1, 1/q1 ]);
        okc := true; C0 := 0;
        try C0 := HyperellipticCurveFromIgusaClebsch(ICq); catch e okc := false; end try;
        if not okc then
            if Verbose then printf "GLUE: HyperellipticCurveFromIgusaClebsch failed\n"; end if;
            continue;
        end if;
        // size control WITHOUT any factoring (PARI hyperellred); then twist fix
        try C0 := hyperellred(C0); catch e; end try;
        f0, h0 := HyperellipticPolynomials(C0);
        g0 := 4*f0 + h0^2;
        good := 0; best := 0; found := false;
        for d in tws do
            okd := true; nchk := 0;
            for p in Keys(aE) do
                if p le 3 or (ZZ!Discriminant(g0)*d) mod p eq 0 then continue; end if;
                gp := PolynomialRing(GF(p))!(d*g0);
                np := #HyperellipticCurve(gp) - 1 - p;   // -a_p
                if np ne -(aE[p]+aF[p]) then okd := false; break; end if;
                nchk +:= 1;
            end for;
            if okd and nchk ge 5 then found := true; best := d; break; end if;
        end for;
        if not found then
            if Verbose then printf "GLUE: no twist matched (M=[%o,%o;%o,%o])\n", m11,m12,m21,m22; end if;
            continue;
        end if;
        Cfin := QuadraticTwist(C0, best);
        try Cfin := hyperellred(Cfin); catch e; end try;
        gk := Sprintf("%o", G2Invariants(Cfin));
        if gk in seen then continue; end if;
        Include(~seen, gk);
        Append(~out, <Cfin, [m11,m12,m21,m22], best>);
        if Verbose then printf "GLUECURVE twist=%o C: %o\n", best, HyperellipticPolynomials(Cfin); end if;
    end for;
    if Verbose then printf "GLUESTATS N=%o: %o invertible M, %o Lagrangian, %o rational\n", N, nM, nlag, nrat; end if;
    return out;
end function;
