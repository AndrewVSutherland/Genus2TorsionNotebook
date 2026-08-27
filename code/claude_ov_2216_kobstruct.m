//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_kobstruct.m   -- lane OV/2216
//
//  WHY do the rational points of the norm surfaces fail to halve?
//
//  The halving criterion for D = Q+T on the odd M_1(8,2,2) model is
//      lambda0 * delta'_j  square in Q  (j = 1,2,3)      [the SURFACE]
//      lambda0 * delta'_K  square in K                   [the OBSTRUCTION]
//  and the C sweeps only impose the surface conditions plus the purely
//  rational NECESSARY shadow  N(lambda0 delta'_K) = square.
//
//  For K = Q(sqrt(DK)) and xi = alpha + beta*sqrt(DK) in K* with
//  N(xi) = alpha^2 - DK*beta^2 = c^2 a rational square, one has
//      xi in (K*)^2  <=>  (alpha+c)/2  or  (alpha-c)/2  is in (Q*)^2,
//  and  ((alpha+c)/2)*((alpha-c)/2) = DK*(beta/2)^2, so the two square
//  classes are  eps  and  DK/eps.  Hence a SINGLE invariant decides it:
//
//      eps(P,T) := squarefree part of (alpha+c)/2  in Q*/(Q*)^2,
//      D = Q+T is 2-divisible  <=>  eps = 1  or  eps = DK   (mod squares).
//
//  This script prints eps for every candidate point and every one of the
//  eight twists, and cross-validates the delta verdict against Magma's
//  own IsDivisibleBy on the Jacobian.
//
//  Usage:
//    code/claude_magma_slot.sh -b candidate_file:=data/xxx.txt \
//        code/claude_ov_2216_kobstruct.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "code/claude_ov_2216_delta.m";

if not assigned candidate_file then
    candidate_file := "data/claude_ov_2216_s1_H2000_uv.txt";
end if;
if not assigned MemGB then MemGB := 6;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Qq := Rationals();
P<X> := PP;

function ParseQ(s)
    if "/" in s then
        pr := Split(s, "/");
        return Qq!StringToInteger(pr[1]) / Qq!StringToInteger(pr[2]);
    end if;
    return Qq!StringToInteger(s);
end function;

lines := Split(Read(candidate_file), "\n");
printf "claude_ov_2216_kobstruct  candidate_file=%o\n", candidate_file;
printf "# eps = squarefree((alpha+c)/2);  2-divisible <=> eps in {1, DK}\n";

npt := 0; ndiv := 0; nsurf := 0;
epsseen := {};

for raw in lines do
    s := raw;
    if #s gt 0 and s[#s] eq "\r" then s := s[1..#s-1]; end if;
    if #s eq 0 or s[1] eq "#" then continue; end if;
    prt := [t : t in Split(s, " ") | #t gt 0];
    if #prt lt 2 then continue; end if;
    u := ParseQ(prt[1]);  v := ParseQ(prt[2]);

    Bc := u^2*v - u^2 + u*v^2 - u - v^2 - v - 2;
    Cc := u^2 + u*v + v^2 + u + v + 1;
    qtm := X^2 - Bc*X + Cc;
    f := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*(-qtm);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        printf "SKIP singular u=%o v=%o\n", u, v; continue;
    end if;
    if not IsIrreducible(qtm) then
        printf "SKIP qtm reducible u=%o v=%o\n", u, v; continue;
    end if;
    L := 1;
    for i in [0..5] do L := LCM(L, Denominator(Coefficient(f, i))); end for;
    fI := P!(L^2*f);
    roots := [ 1/(u-1), 1/(v-1), -1/(u+v+2) ];
    for r in roots do assert Evaluate(fI, r) eq 0; end for;

    KK<th> := NumberField(qtm);
    DK := Bc^2 - 4*Cc;                    // sqrt(DK) = 2*th - Bc
    dksq := CO_SqClass(DK);

    C2 := HyperellipticCurve(fI);
    J := Jacobian(C2);
    yq := u*v*(u+v+1);
    assert (L*yq)^2 eq Evaluate(fI, Qq!-1);
    DQ := J![X + 1, Qq!(L*yq)];
    ok, dpQ := CO_DeltaPrimeCoprime(X+1, roots, KK, th);
    assert ok;

    lin := [ X - roots[i] : i in [1..3] ];
    uSlist := [ P!1, lin[1]*lin[2], lin[1]*lin[3], lin[2]*lin[3], qtm,
                lin[2]*lin[3]*qtm, lin[1]*lin[3]*qtm, lin[1]*lin[2]*qtm ];
    tlab   := [ "0     ", "A_r1r2", "A_r1r3", "A_r2r3", "B_quad",
                "C_r1  ", "C_r2  ", "C_r3  " ];
    JA12 := J![lin[1]*lin[2], 0];
    JA13 := J![lin[1]*lin[3], 0];
    JA23 := J![lin[2]*lin[3], 0];
    JBq  := J![qtm, 0];
    JT := [ J!0, JA12, JA13, JA23, JBq, JA23+JBq, JA13+JBq, JA12+JBq ];

    npt +:= 1;
    printf "\nPOINT u=%o v=%o  DK_sqclass=%o  L=%o\n", u, v, dksq, L;

    for k in [1..8] do
        if k eq 1 then
            dp := dpQ;
        else
            ok2, dpT := CO_DeltaPrimeTwo(uSlist[k], fI, roots, KK, th);
            if not ok2 then printf "  %o  DELTA_FAIL\n", tlab[k]; continue; end if;
            dp := CO_MulPrime(dpQ, dpT);
        end if;
        nrm := dp[1]*dp[2]*dp[3]*Norm(dp[4]);
        lam := CO_SqClass(nrm);
        rat1 := IsSquare(lam*dp[1]);
        rat2 := IsSquare(lam*dp[2]);
        rat3 := IsSquare(lam*dp[3]);
        ratOK := rat1 and rat2 and rat3;

        xi := (KK!lam)*dp[4];
        // xi = a + b*th = alpha + beta*sqrt(DK),  alpha = a + b*Bc/2, beta = b/2
        cf := Eltseq(xi);
        aa := Qq!cf[1]; bb := Qq!cf[2];
        alpha := aa + bb*Bc/2;
        beta  := bb/2;
        Nxi := Norm(xi);
        issq, c := IsSquare(Nxi);
        epsp := 0; epsm := 0;
        if issq then
            epsp := CO_SqClass((alpha + c)/2 eq 0 select 1 else (alpha + c)/2);
            epsm := CO_SqClass((alpha - c)/2 eq 0 select 1 else (alpha - c)/2);
        end if;
        ksq := IsSquare(xi);
        DT := DQ + JT[k];
        ordDT := Order(DT);
        magdiv := IsDivisibleBy(DT, 2);
        deltadiv := ratOK and ksq;
        assert deltadiv eq magdiv;          // delta module vs Magma, per twist

        if ratOK then
            nsurf +:= 1;
            Include(~epsseen, epsp);
        end if;
        if magdiv then ndiv +:= 1; end if;

        printf "  T=%o ord=%o  ratsq=(%o,%o,%o) Nxi_sq=%o eps+=%o eps-=%o  DK=%o  Ksq=%o  divisible=%o\n",
            tlab[k], ordDT, rat1, rat2, rat3, issq, epsp, epsm, dksq, ksq, magdiv;
    end for;
end for;

printf "\npoints=%o  (point,twist) pairs on a norm surface=%o  2-divisible=%o\n",
       npt, nsurf, ndiv;
printf "EPS_CLASSES_SEEN_ON_SURFACE %o\n", Sort(SetToSequence(epsseen));
print "SEARCH_DONE";
quit;
