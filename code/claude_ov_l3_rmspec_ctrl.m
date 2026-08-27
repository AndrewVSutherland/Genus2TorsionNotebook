// ===========================================================================
// Lane 3 : the CONDITIONED control for claude_ov_l3_rmspectrum.m.
//
// The RM-spectrum test histograms the squarefree core of  n = a1^2 - 4(a2-2p)
// over Sigma''(F_p) and looks for a spike (= a component of Sigma'' with RM).
// But every point of Sigma'' has  11 | #J(F_p) = L(1),  which by itself biases
// the pair (a1,a2) and hence the distribution of n.  So the raw ambient
// histogram is the WRONG control; the right one is
//
//     random points of the chart CONDITIONED on 11 | #J(F_p).
//
// This script produces it.  Output lines "B <core>" (11 | #J) and "N <core>"
// (11 does not divide #J), from which both the conditioned and unconditioned
// ambient histograms can be read.
//
// usage: code/claude_magma_slot.sh -b P:=101 NTRY:=24000 NCH:=8 code/claude_ov_l3_rmspec_ctrl.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned P then P := 101; elif Type(P) eq MonStgElt then P := StringToInteger(P); end if;
if not assigned NTRY then NTRY := 24000; elif Type(NTRY) eq MonStgElt then NTRY := StringToInteger(NTRY); end if;
if not assigned NCH then NCH := 8; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned OUT then OUT := Sprintf("results/claude_ov_l3_rmspec_ctrl_p%o.log", P); end if;

Z := Integers();
Fp := GF(P);
Pp<x> := PolynomialRing(Fp);

dir := "/tmp/claude_ov_l3_rmctrl";
System("mkdir -p " cat dir cat " && rm -f " cat dir cat "/part*.txt");
for ch in [0..NCH-1] do
    pid := Fork();
    if pid eq 0 then
        SetSeed(70000+ch);
        F := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
        cnt := 0;
        while cnt lt (NTRY div NCH) do
            a := Random(Fp); b := Random(Fp); c := Random(Fp); d := Random(Fp);
            f := x*(x-1)*(x^2+a*x+b)*(x^2+c*x+d);
            if Degree(f) ne 6 or not IsSquarefree(f) then continue; end if;
            C := HyperellipticCurve(f);
            L := LPolynomial(C);
            a1 := Z!(-Coefficient(L,1));  a2 := Z!Coefficient(L,2);
            nJ := 1 - a1 + a2 - P*a1 + P^2;
            n := a1^2 - 4*(a2 - 2*P);
            core := (n eq 0) select 0 else SquarefreeFactorization(n);
            Puts(F, Sprintf("%o %o", (nJ mod 11 eq 0) select "B" else "N", core));
            cnt +:= 1;
        end while;
        Flush(F); delete F;
        quit;
    end if;
end for;
WaitForAllChildren();

hB := AssociativeArray(); hN := AssociativeArray(); nB := 0; nN := 0;
for ch in [0..NCH-1] do
    for l in Split(Read(dir cat "/part" cat IntegerToString(ch) cat ".txt"), "\n") do
        if #l eq 0 then continue; end if;
        s := Split(l, " ");  core := StringToInteger(s[2]);
        if s[1] eq "B" then
            nB +:= 1; if not IsDefined(hB, core) then hB[core] := 0; end if; hB[core] +:= 1;
        else
            nN +:= 1; if not IsDefined(hN, core) then hN[core] := 0; end if; hN[core] +:= 1;
        end if;
    end for;
end for;

G := Open(OUT, "w");
procedure Say(G, s) printf "%o\n", s; Puts(G, s); end procedure;
Say(G, Sprintf("RMSPEC_CTRL p=%o  with 11|#J : %o   without : %o", P, nB, nN));
procedure Top(G, h, n, tag)
    ks := Sort(Setseq(Keys(h)), func<u,v | h[v]-h[u]>);
    Say(G, Sprintf("-- %o (n=%o), top 12:", tag, n));
    for i in [1..Minimum(12,#ks)] do
        Say(G, Sprintf("   core=%-6o count=%-6o frac=%o", ks[i], h[ks[i]],
                       RealField(4)!(h[ks[i]]/n)));
    end for;
end procedure;
Top(G, hB, nB, "AMBIENT | 11 divides #J(F_p)   <-- the correct control");
Top(G, hN, nN, "AMBIENT | 11 does not divide #J(F_p)");
Flush(G); delete G;
printf "RMCTRL_DONE\n";
quit;
