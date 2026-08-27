//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_ktest.m
//
//  Downstream stage for code/claude_ov_2216_sweep.c.
//
//  Reads lines  "HIT <surface> <p> <q> <a> <b>"  (u = p/q, k = a/b),
//  reconstructs v from the conic parametrization, then
//    1. EXACTLY re-tests the rational square condition (the C stage only
//       applies a quadratic-residue prefilter, so a few percent of the
//       input can be false positives),
//    2. applies the quadratic-field ("K") condition with the validated
//       x-T module, which decides  [2,2,16] <= J(Q)  outright,
//    3. runs exact TorsionSubgroup on anything that survives.
//
//  Usage:
//    magma -b infile:=results/xxx.txt code/claude_ov_2216_ktest.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "code/claude_ov_2216_delta.m";

if not assigned infile then infile := "results/claude_ov_2216_sweep_s0.txt"; end if;
if not assigned MemGB then MemGB := 6;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

P<X> := PP;

// reconstruct v from (surface, u, k)
//   surface 0:  v =  A/(A - k^2)
//   surface 1:  v = -A(2u+1)/(k^2 + A)     with A = u(u-1)
function VFrom(surface, uu, kk)
    A := uu*(uu-1);
    if surface eq 0 then
        den := A - kk^2;
        if den eq 0 then return false, _; end if;
        return true, A/den;
    else
        den := kk^2 + A;
        if den eq 0 then return false, _; end if;
        return true, -A*(2*uu+1)/den;
    end if;
end function;

lines := Split(Read(infile), "\n");
printf "claude_ov_2216_ktest infile=%o lines=%o\n", infile, #lines;

nrow := 0; nrec := 0; nondeg := 0; irred := 0;
nY := 0; nfull := 0; nfalse := 0;
hits := [];

for raw in lines do
    s := raw;
    if #s gt 0 and s[#s] eq "\r" then s := s[1..#s-1]; end if;
    if #s eq 0 or s[1] eq "#" then continue; end if;
    parts := [t : t in Split(s, " ") | #t gt 0];
    if #parts lt 5 or parts[1] ne "HIT" then continue; end if;
    nrow +:= 1;
    surface := StringToInteger(parts[2]);
    p := StringToInteger(parts[3]); q := StringToInteger(parts[4]);
    a := StringToInteger(parts[5]); b := StringToInteger(parts[6]);
    uu := Qq!p/q; kk := Qq!a/b;
    okv, vv := VFrom(surface, uu, kk);
    if not okv then continue; end if;
    nrec +:= 1;
    u := uu; v := vv;
    if u eq 0 or v eq 0 or u eq v or u eq 1 or v eq 1 then continue; end if;
    if u+v+1 eq 0 or u+v+2 eq 0 then continue; end if;

    B := u^2*v - u^2 + u*v^2 - u - v^2 - v - 2;
    C := u^2 + u*v + v^2 + u + v + 1;
    qt := -X^2 + B*X - C;
    f0 := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qt;
    if Degree(f0) ne 5 or Discriminant(f0) eq 0 then continue; end if;
    nondeg +:= 1;
    // q irreducible  <=>  factor type [1,1,1,2]  <=>  2-rank 3  <=> EXACT
    // [2,2,16] is possible; if it splits the target would be [2,2,2,16].
    qsplit := not IsIrreducible(qt);
    if not qsplit then irred +:= 1; end if;
    if qsplit then
        printf "QSPLIT u=%o v=%o (would give [2,2,2,16], order 128)\n", u, v;
        continue;
    end if;

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

    dl := [* dQ *];
    for k in [1..3] do
        uS := qtm * &*[ lin[j] : j in [1..3] | j ne k ];
        okk, dT := CO_DeltaPrimeTwo(uS, f, roots, KK, th);
        if not okk then continue; end if;
        Append(~dl, CO_MulPrime(dQ, dT));
    end for;
    if #dl ne 4 then continue; end if;

    anyY := false;
    for idx in [1..4] do
        dp := dl[idx];
        s1 := CO_SqClass(dp[1]);
        if CO_SqClass(dp[2]) ne s1 then continue; end if;
        if CO_SqClass(dp[3]) ne s1 then continue; end if;
        anyY := true; nY +:= 1;
        printf "YPOINT surface=%o class=%o u=%o v=%o d=%o\n", surface, idx-1, u, v, s1;
        if not IsSquare((KK!s1)*dp[4]) then continue; end if;
        nfull +:= 1;
        printf "FULLTRIVIAL surface=%o class=%o u=%o v=%o d=%o\n", surface, idx-1, u, v, s1;
        // exact confirmation
        Cc := HyperellipticCurve(f);
        J := Jacobian(Cc);
        G, phi := TorsionSubgroup(J);
        printf "  EXACT TORSION u=%o v=%o -> %o\n", u, v, Invariants(G);
        Append(~hits, <u, v, idx-1, Invariants(G)>);
    end for;
    if not anyY then nfalse +:= 1; end if;
end for;

printf "rows=%o reconstructed=%o nondegenerate=%o qt_irreducible=%o\n",
       nrow, nrec, nondeg, irred;
printf "Ypoints=%o fulltrivial=%o prefilter_false_positives=%o\n", nY, nfull, nfalse;
printf "hits=%o\n", #hits;
for H in hits do print "HIT", H; end for;
print "SEARCH_DONE";
quit;
