// glue_window.m — construct the (N,N)-gluing of an m-isogenous pair in a
// golden window ((N,m) in {(5,6),(5,24),(7,10)}): psi = c*phi|E[N] is
// Galois-equivariant (isogeny-induced) and Kani-NONdegenerate (p-adic
// descents; d^2 = N(a^2+m b^2) has no solutions), so A = (ExF)/graph(psi)
// IS a Jacobian over Q and T1 x T2 injects into its rational torsion.
// The matrix of phi on H_1 is computed numerically from the isogeny's
// period scalar (phi corresponds to z -> alpha z with alpha*Lambda_E
// inside Lambda_F), giving the single correct M; tau -> Rosenhain ->
// invariants (rational by theory) -> Mestre -> twist -> certification.
//
// Inputs (assigned or defaults): COND, CLSIDX, K1, K2, NN, prec
// Usage example:
//   magma -b COND:=302 CLSIDX:=1 K1:=1 K2:=2 NN:=5 glue_window.m
SetColumns(0);
SetMemoryLimit(8*10^9);
if not assigned COND then error "need COND"; elif Type(COND) eq MonStgElt then COND := StringToInteger(COND); end if;
if not assigned CLSIDX then CLSIDX := 1; elif Type(CLSIDX) eq MonStgElt then CLSIDX := StringToInteger(CLSIDX); end if;
if not assigned K1 then error "need K1"; elif Type(K1) eq MonStgElt then K1 := StringToInteger(K1); end if;
if not assigned K2 then error "need K2"; elif Type(K2) eq MonStgElt then K2 := StringToInteger(K2); end if;
if not assigned NN then NN := 5; elif Type(NN) eq MonStgElt then NN := StringToInteger(NN); end if;
if not assigned prec then prec := 300; elif Type(prec) eq MonStgElt then prec := StringToInteger(prec); end if;
load "split_lab.m";
load "analytic_glue.m";

QQ := Rationals(); ZZ := Integers();
D := CremonaDatabase();
cls := EllipticCurves(D, COND, CLSIDX);
EA := cls[K1]; EB := cls[K2];
T1 := Invariants(TorsionSubgroup(EA)); T2 := Invariants(TorsionSubgroup(EB));
okI, phi := IsIsogenous(EA, EB);
error if not okI, "not isogenous";
m := Degree(phi);
printf "pair %o x %o, torsions %o %o, deg(phi)=%o, N=%o\n",
    aInvariants(EA), aInvariants(EB), T1, T2, m, NN;
// window check: c^2*m = -1 mod N solvable
cs := [ c : c in [1..NN-1] | (c^2*m + 1) mod NN eq 0 ];
error if IsEmpty(cs), "no anti-isometric scaling: not in a window";
printf "anti scalings c = %o\n", cs;

CC := ComplexField(prec);
wE := NormPeriods(EA, CC); wF := NormPeriods(EB, CC);
// find alpha with alpha*Lambda_E subset Lambda_F of index m:
// solve alpha*wE[i] = R[i][1]*wF[1] + R[i][2]*wF[2] with integers R, det(R) = +-m.
// alpha is determined up to Aut = +-1; get it from the two lattice equations:
// try candidate alpha = (p*wF[1]+q*wF[2])/wE[1] over small integers p,q and
// test that alpha*wE[2] is also in Lambda_F.
function InLat(z, w1, w2, eps)
    // solve z = x*w1 + y*w2 over R^2, check near-integrality
    M := Matrix(RealField(Precision(Parent(z))), 2,2,
        [ Re(w1), Re(w2), Im(w1), Im(w2) ]);
    v := Vector(RealField(Precision(Parent(z))), [ Re(z), Im(z) ]);
    s := Solution(Transpose(M), v);
    x := s[1]; y := s[2];
    xr := Round(x); yr := Round(y);
    if Abs(x - xr) lt eps and Abs(y - yr) lt eps then return true, xr, yr; end if;
    return false, 0, 0;
end function;
eps := RealField(20)!10.0^(-prec div 3);
found := false; alpha := CC!0; R := Matrix(ZZ,2,2,[0,0,0,0]);
for p in [-m..m] do
    for q in [-m..m] do
        if p eq 0 and q eq 0 then continue; end if;
        al := (p*wF[1] + q*wF[2])/wE[1];
        ok2, r21, r22 := InLat(al*wE[2], wF[1], wF[2], eps);
        if not ok2 then continue; end if;
        dt := p*r22 - q*r21;
        if Abs(dt) ne m then continue; end if;
        found := true; alpha := al;
        R := Matrix(ZZ,2,2,[p,q,r21,r22]);
        break p;
    end for;
end for;
error if not found, "period scalar for phi not found (widen search)";
printf "phi matrix on H1 (rows = images of wE-basis in wF-basis): %o det %o\n", R, Determinant(R);

for c in cs do
    // psi = c*phi on N-torsion: x = (a*wE1+b*wE2)/N maps to c*(a,b)*R * (wF/N)
    MR := Matrix(Integers(NN), 2, 2, [ c*R[1][1], c*R[1][2], c*R[2][1], c*R[2][2] ]);
    m11 := ZZ!MR[1][1]; m12 := ZZ!MR[1][2]; m21 := ZZ!MR[2][1]; m22 := ZZ!MR[2][2];
    printf "== c=%o: M = [%o,%o;%o,%o] (det mod N = %o) ==\n", c, m11,m12,m21,m22,
        (m11*m22-m12*m21) mod NN;
    rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
              [QQ|1/NN, 0, m11/NN, m12/NN], [QQ|0, 1/NN, m21/NN, m22/NN] ];
    H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| NN*x : x in r ] : r in rows ]));
    B := [ [ QQ | H[i][j]/NN : j in [1..4] ] : i in [1..4] ];
    J := Matrix(QQ, 4,4, [ [ NN*EProdPair(B[i], B[j]) : j in [1..4] ] : i in [1..4] ]);
    if not forall{ <i,j> : i,j in [1..4] | Denominator(J[i][j]) eq 1 } then
        printf "   not Lagrangian (unexpected) -- skip\n"; continue;
    end if;
    JZ := Matrix(ZZ, 4,4, [ [ ZZ!J[i][j] : j in [1..4] ] : i in [1..4] ]);
    F0, T := FrobeniusFormAlternating(JZ);
    Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
    cols := [ [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ] where v := Cb[i] : i in [1..4] ];
    PA := Matrix(CC, 2,2, [ cols[1][1], cols[2][1], cols[1][2], cols[2][2] ]);
    PB := Matrix(CC, 2,2, [ cols[3][1], cols[4][1], cols[3][2], cols[4][2] ]);
    tau := PA^-1*PB;
    if Im(tau[1][1]) lt 0 then tau := PB^-1*PA; end if;
    okr := true; ros := [];
    try ros := RosenhainInvariants(tau); catch e okr := false; end try;
    if not okr then printf "   rosenhain failed\n"; continue; end if;
    PCx<xx> := PolynomialRing(CC);
    g := xx*(xx-1)*&*[ xx - r : r in ros ];
    IC := IgusaClebschInvariants(HyperellipticCurve(g));
    if Abs(IC[4]) lt 10.0^(-40) then printf "   DEGENERATE (I10 ~ 0)?!\n"; continue; end if;
    j1 := IC[1]^5/IC[4]; j2 := IC[1]^3*IC[2]/IC[4]; j3 := IC[1]^2*IC[3]/IC[4];
    qs := [ QQ | ]; okall := true;
    for jj in [j1,j2,j3] do
        r := Re(jj);
        okq := false;
        for hb in [20, 40, 60, 90] do
            if 2*hb + 20 gt prec then break; end if;
            q := BestApproximation(r, 10^hb);
            hq := Max(Ilog(10, 1+Abs(Numerator(q))), Ilog(10, 1+Denominator(q)));
            if Abs(r - q) lt RealField(20)!10.0^(-(2*hq + 15)) then
                Append(~qs, q); okq := true; break;
            end if;
        end for;
        if not okq then okall := false; break; end if;
    end for;
    if not okall then
        printf "   invariants NOT verified rational at this precision (heights too large?)\n";
        printf "   j1 approx %o\n", ComplexField(30)!j1;
        continue;
    end if;
    printf "   RATIONAL INVARIANTS verified: heights %o digits\n",
        [ Ilog(10,1+Abs(Numerator(q))) + Ilog(10,1+Denominator(q)) : q in qs ];
    ICq := ReduceIC([ QQ | 1, qs[2]/qs[1], qs[3]/qs[1], 1/qs[1] ]);
    printf "   reduced IC heights %o\n", [ Ilog(10, 1+Abs(Numerator(x))) + Ilog(10, 1+Denominator(x)) : x in ICq ];
    t0 := Cputime();
    okc := true; C0 := 0;
    try C0 := HyperellipticCurveFromIgusaClebsch(ICq); catch e okc := false; end try;
    if not okc then printf "   Mestre failed\n"; continue; end if;
    printf "   MESTRE OK %o s\n", Cputime()-t0;
    try C0 := hyperellred(C0); catch e; end try;
    f0, h0 := HyperellipticPolynomials(C0);
    g0 := 4*f0 + h0^2;
    tws := [ ZZ | d : d in Divisors(SquarefreeFactorization(2*COND*NN)) ];
    tws := tws cat [-d : d in tws];
    foundtw := false; best := 0;
    for d in tws do
        okd := true; nchk := 0;
        for p in PrimesUpTo(80) do
            if (2*ZZ!Discriminant(g0)*d*COND*NN) mod p eq 0 then continue; end if;
            gp := PolynomialRing(GF(p))!(d*g0);
            np := #HyperellipticCurve(gp) - 1 - p;
            if np ne -(TraceOfFrobenius(EA,p)+TraceOfFrobenius(EB,p)) then okd := false; break; end if;
            nchk +:= 1;
        end for;
        if okd and nchk ge 5 then foundtw := true; best := d; break; end if;
    end for;
    if not foundtw then printf "   no twist matched\n"; continue; end if;
    C1 := QuadraticTwist(C0, best);
    try C1 := hyperellred(C1); catch e; end try;
    okS, gg := IntegralSextic(C1);
    if not okS then printf "   bad model\n"; continue; end if;
    einv := [ZZ|-1];
    try einv := ExactTorsion(gg); catch e printf "   torsion failed\n"; continue; end try;
    prod := Invariants(AbelianGroup(T1 cat T2));
    printf "WINDOW_RESULT cond=%o c=%o twist=%o invs=%o (contains %o expected) %o g=%o\n",
        COND, c, best, einv, prod, (not einv in KNOWN) select "*** NEW GROUP ***" else "known", gg;
end for;
printf "GLUE_WINDOW_DONE\n";
quit;
