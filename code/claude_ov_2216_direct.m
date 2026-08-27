//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_direct.m
//
//  GROUND-TRUTH direct enumeration over (u,v) of bounded height on the
//  odd M_1(8,2,2) chart.  For every nondegenerate (u,v) it computes the
//  delta' vectors of the four genuinely distinct order-8 classes
//        Q,  Q+T_1,  Q+T_2,  Q+T_3
//  with the VALIDATED module code/claude_ov_2216_delta.m, and records
//    (a) "Y-points": the three RATIONAL components share a common square
//        class  (i.e. the point lies on the norm surface Y_t), and
//    (b) full delta-triviality (Y-point AND the quadratic-field
//        condition) -- which is equivalent to [2,2,16] <= J(Q).
//
//  This is the reference against which the C sweeper
//  code/claude_ov_2216_sweep.c is checked for completeness, and it
//  independently reproduces the recorded negatives.
//
//  Usage: magma -b height:=H code/claude_ov_2216_direct.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "code/claude_ov_2216_delta.m";

if not assigned height then height := 30;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned MemGB then MemGB := 4;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

P<X> := PP;

hs := [];
for a in [-height..height] do
    for b in [1..height] do
        if GCD(AbsoluteValue(a), b) eq 1 then Append(~hs, Qq!a/b); end if;
    end for;
end for;
hs := Sort(Setseq(Seqset(hs)));

printf "claude_ov_2216_direct height=%o  #rationals=%o\n", height, #hs;

checked := 0; nondeg := 0; irred := 0;
ypts := [0,0,0,0];      // Y-points per class (untwisted, T1, T2, T3)
full := [0,0,0,0];
hits := [];

for u in hs do
    for v in hs do
        checked +:= 1;
        if u eq 0 or v eq 0 or u eq v or u eq 1 or v eq 1 then continue; end if;
        if u+v+1 eq 0 or u+v+2 eq 0 then continue; end if;
        B := u^2*v - u^2 + u*v^2 - u - v^2 - v - 2;
        C := u^2 + u*v + v^2 + u + v + 1;
        qt := -X^2 + B*X - C;
        f0 := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qt;
        if Degree(f0) ne 5 then continue; end if;
        if Discriminant(f0) eq 0 then continue; end if;
        nondeg +:= 1;
        if not IsIrreducible(qt) then continue; end if;   // want type [1,1,1,2]
        irred +:= 1;

        Lden := 1;
        for i in [0..5] do Lden := LCM(Lden, Denominator(Coefficient(f0,i))); end for;
        f := P!(Lden^2*f0);
        roots := [ Qq!rt[1] : rt in Roots(f) ];
        if #roots ne 3 then continue; end if;
        qtm := qt/LeadingCoefficient(qt);
        KK<th> := NumberField(qtm);
        lin := [ X - r : r in roots ];

        ok, dQ := CO_DeltaPrimeCoprime(X+1, roots, KK, th);
        if not ok then continue; end if;

        // the three twists T_k = [(r_k,0)-inf], encoded by the
        // complementary even finite subset  S = {r_i,r_j,th,th'}
        dl := [* dQ *];
        for k in [1..3] do
            uS := qtm * &*[ lin[j] : j in [1..3] | j ne k ];
            okk, dT := CO_DeltaPrimeTwo(uS, f, roots, KK, th);
            if not okk then continue; end if;
            Append(~dl, CO_MulPrime(dQ, dT));
        end for;
        if #dl ne 4 then continue; end if;

        for idx in [1..4] do
            dp := dl[idx];
            s1 := CO_SqClass(dp[1]);
            if CO_SqClass(dp[2]) ne s1 then continue; end if;
            if CO_SqClass(dp[3]) ne s1 then continue; end if;
            ypts[idx] +:= 1;
            printf "YPOINT class=%o u=%o v=%o d=%o\n", idx-1, u, v, s1;
            if not IsSquare((KK!s1)*dp[4]) then continue; end if;
            full[idx] +:= 1;
            printf "FULLTRIVIAL class=%o u=%o v=%o d=%o\n", idx-1, u, v, s1;
            Append(~hits, <u, v, idx-1>);
        end for;
    end for;
end for;

printf "checked=%o nondegenerate=%o qt_irreducible=%o\n", checked, nondeg, irred;
printf "YPOINTS untwisted=%o T1=%o T2=%o T3=%o\n", ypts[1], ypts[2], ypts[3], ypts[4];
printf "FULL    untwisted=%o T1=%o T2=%o T3=%o\n", full[1], full[2], full[3], full[4];
printf "hits=%o\n", #hits;
for H in hits do print "HIT", H; end for;
print "SEARCH_DONE";
quit;
