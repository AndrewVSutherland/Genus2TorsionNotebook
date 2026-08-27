// ===========================================================================
// Lane 3 : is Sigma'' contained in SOME Humbert surface other than H_5?
//
// Refuting containment in H_5 (claude_ov_l3_humbtest.m) leaves one loophole:
// Sigma'' could still be RM-forced by a DIFFERENT discriminant (another real
// quadratic field, or a non-maximal order of Q(sqrt 5)).  Deriving H_Delta for
// each Delta needs an Elkies-Kumar parametrisation we do not have, so test it
// intrinsically instead:
//
//   for a point of Sigma''(F_p), the real Weil subfield of Q(pi) is
//   Q(sqrt n),  n = a1^2 - 4(a2 - 2p).   If a whole COMPONENT of Sigma'' had
//   RM by discriminant Delta, then at every F_p-point of that component (where
//   the reduction is ordinary and F_p-simple and the RM is F_p-rational) the
//   squarefree core of n would be the core of Delta.  A component is ~1/8.3 of
//   Sigma''(F_p), so an RM component shows up as a ~12% SPIKE in the histogram
//   of cores -- against an ambient background in which each core occurs O(1/p).
//
// Built-in positive control: the points of Sigma''(F_p) that lie ON H_5 must
// spike at core = 5.
//
// usage: code/claude_magma_slot.sh -b P:=101 NS:=9000 NCH:=4 code/claude_ov_l3_rmspectrum.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned P then P := 101; elif Type(P) eq MonStgElt then P := StringToInteger(P); end if;
if not assigned NS then NS := 9000; elif Type(NS) eq MonStgElt then NS := StringToInteger(NS); end if;
if not assigned NCH then NCH := 4; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned DUMP then DUMP := Sprintf("/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/sigma_p%o.txt", P); end if;
if not assigned OUT then OUT := Sprintf("results/claude_ov_l3_rmspectrum_p%o.log", P); end if;

Z := Integers();  Q := Rationals();
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5Q := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");
Fp := GF(P);
Rp := PolynomialRing(Fp, 4);
H5p := Rp ! H5Q;
Pp<x> := PolynomialRing(Fp);

// returns  <ok, onH5, core>   ok = 0 if degenerate
function Probe(a,b,c,d)
    f := x*(x-1)*(x^2+a*x+b)*(x^2+c*x+d);
    if Degree(f) ne 6 or not IsSquarefree(f) then return 0, 0, 0; end if;
    C := HyperellipticCurve(f);
    I := IgusaClebschInvariants(C);
    on := (Evaluate(H5p, [I[1],I[2],I[3],I[4]]) eq 0) select 1 else 0;
    L := LPolynomial(C);
    a1 := Z!(-Coefficient(L,1));  a2 := Z!Coefficient(L,2);
    n := a1^2 - 4*(a2 - 2*P);
    if n eq 0 then return 1, on, 0; end if;
    return 1, on, SquarefreeFactorization(n);
end function;

// AGGONLY:=1 re-does only the aggregation, reusing the existing part files
if not assigned AGGONLY then AGGONLY := 0;
elif Type(AGGONLY) eq MonStgElt then AGGONLY := StringToInteger(AGGONLY); end if;

dir := "/tmp/claude_ov_l3_rmspec";
if AGGONLY eq 0 then
lines := [l : l in Split(Read(DUMP), "\n") | #l gt 2];
printf "SIGMA_POINTS %o  sampling %o\n", #lines, NS;
step := Maximum(1, #lines div NS);

System("mkdir -p " cat dir cat " && rm -f " cat dir cat "/part*.txt");
for ch in [0..NCH-1] do
    pid := Fork();
    if pid eq 0 then
        F := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
        // ---- Sigma'' sample
        idx := 1 + ch*step;
        while idx le #lines do
            s := Split(lines[idx], " ");
            ok, on, core := Probe(Fp!StringToInteger(s[1]), Fp!StringToInteger(s[2]),
                                  Fp!StringToInteger(s[3]), Fp!StringToInteger(s[4]));
            if ok eq 1 then Puts(F, Sprintf("S %o %o", on, core)); end if;
            idx +:= NCH*step;
        end while;
        // ---- ambient control, same size
        SetSeed(1000+ch);
        cnt := 0;
        while cnt lt (NS div NCH) do
            ok, on, core := Probe(Random(Fp),Random(Fp),Random(Fp),Random(Fp));
            if ok eq 1 then Puts(F, Sprintf("A %o %o", on, core)); cnt +:= 1; end if;
        end while;
        Flush(F); delete F;
        quit;
    end if;
end for;
WaitForAllChildren();
end if;

hS := AssociativeArray(); hSoff := AssociativeArray();
hSon := AssociativeArray(); hA := AssociativeArray();
nS := 0; nSon := 0; nSoff := 0; nA := 0;
for ch in [0..NCH-1] do
    for l in Split(Read(dir cat "/part" cat IntegerToString(ch) cat ".txt"), "\n") do
        if #l gt 2 then
            s := Split(l, " ");
            on := StringToInteger(s[2]); core := StringToInteger(s[3]);
            if s[1] eq "S" then
                nS +:= 1;
                if not IsDefined(hS, core) then hS[core] := 0; end if;
                hS[core] +:= 1;
                if on eq 1 then
                    nSon +:= 1;
                    if not IsDefined(hSon, core) then hSon[core] := 0; end if;
                    hSon[core] +:= 1;
                else
                    nSoff +:= 1;
                    if not IsDefined(hSoff, core) then hSoff[core] := 0; end if;
                    hSoff[core] +:= 1;
                end if;
            else
                nA +:= 1;
                if not IsDefined(hA, core) then hA[core] := 0; end if;
                hA[core] +:= 1;
            end if;
        end if;
    end for;
end for;

G := Open(OUT, "w");
procedure Say(G, s) printf "%o\n", s; Puts(G, s); end procedure;
Say(G, Sprintf("RMSPECTRUM p=%o  sigma_sampled=%o (onH5 %o, offH5 %o)  ambient=%o",
    P, nS, nSon, nSoff, nA));

procedure Top(G, h, n, tag)
    ks := Sort(Setseq(Keys(h)), func<u,v | h[v]-h[u]>);
    Say(G, Sprintf("-- %o : %o distinct cores over %o points; top 12:", tag, #ks, n));
    for i in [1..Minimum(12,#ks)] do
        Say(G, Sprintf("   core=%-6o count=%-6o frac=%o", ks[i], h[ks[i]],
                       RealField(4)!(h[ks[i]]/n)));
    end for;
end procedure;
Top(G, hS,    nS,    "SIGMA'' (all)");
if nSon  gt 0 then Top(G, hSon,  nSon,  "SIGMA'' cap H_5   <-- positive control"); end if;
if nSoff gt 0 then Top(G, hSoff, nSoff, "SIGMA'' off H_5   <-- the question"); end if;
Top(G, hA,    nA,    "AMBIENT chart (unconditioned)");
Flush(G); delete G;
printf "RMSPEC_DONE\n";
quit;
