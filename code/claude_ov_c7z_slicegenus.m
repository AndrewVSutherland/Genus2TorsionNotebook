// Lane 2 (claude_ov_c7z): which slices z_1 = c of the contact-7 three-root surface are SPECIAL?
// The generic slice is a (3,3) curve of genus 4; c=2 drops to genus 1 (rank 1).
// This forks over a list of rational c and records the genus of every irreducible component.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
Pw<w> := PolynomialRing(Q);

cands := [];
for d in [1..12] do
  for n in [-30..30] do
    if GCD(Abs(n),d) ne 1 then continue; end if;
    cv := n/d;
    if cv in {Q|0,1,1/2} then continue; end if;
    Append(~cands, cv);
  end for;
end for;
// z-values occurring in the 11 known three-root triples (each triple lies on 3 slices)
known := [ Q | -1/2,4,5/2, -1/9,2,-7/3, -1/4,-8/7,22/7, 2,-8/7,19/4, 9/5,17/21,25/21,
           17/21,18/13,49/39, -61/450,625/114,2, -41/124,-16/17,289/124, 297/133,2,361/525,
           50/33,189/155,121/155, 2,49/36,50/63 ];
for k in known do if not k in cands then Append(~cands, k); end if; end for;
printf "SLICEGENUS scanning %o values of c\n", #cands;

// Output directory: pass dir:=... on the command line, else use the
// repository-relative default (run from the repo root).
if not assigned dir then dir := "results/slicegenus"; end if;
System("mkdir -p " cat dir);
NCH := 24;
for ch in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then
    Fh := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
    for i in [ch+1..#cands by NCH] do
      cv := cands[i];
      T := 3 - 1/(1-cv);
      if T eq 1 then Puts(Fh, Sprintf("c=%o DEGENERATE T=1", cv)); continue; end if;
      e1 := 6 - cv; e4 := 2/cv;
      Aa := (T-2)/(T-1);
      Bb := (T*(cv-5+2/cv) - 3*cv + 14)/(T-1);
      R<X,Y> := PolynomialRing(Q,2);
      Nm := func<t| -(t^4 - e1*t^3 - Bb*t + e4)>;
      Dn := func<t| t^2 - Aa*t>;
      N := Nm(X)*Dn(Y) - Nm(Y)*Dn(X);
      okq, FF := IsDivisibleBy(N, X-Y);
      if not okq then Puts(Fh, Sprintf("c=%o NODIV", cv)); continue; end if;
      if FF eq 0 then Puts(Fh, Sprintf("c=%o ZERO", cv)); continue; end if;
      fa := Factorisation(FF);
      gs := [];
      for g in fa do
        gg := -1;
        try gg := Genus(ProjectiveClosure(Curve(AffineSpace(R), g[1]))); catch e gg := -1; end try;
        Append(~gs, <Degree(g[1],X), Degree(g[1],Y), g[2], gg>);
      end for;
      Puts(Fh, Sprintf("c=%o nfac=%o comps=%o", cv, #fa, gs));
    end for;
    Flush(Fh); delete Fh;
    quit;
  end if;
end for;
WaitForAllChildren();
lines := [];
for ch in [0..NCH-1] do
  for l in Split(Read(dir cat "/part" cat IntegerToString(ch) cat ".txt"), "\n") do
    if #l gt 0 then Append(~lines, l); end if;
  end for;
end for;
printf "collected %o lines\n", #lines;
for l in lines do printf "%o\n", l; end for;
printf "SLICEGENUS_DONE n=%o\n", #lines;
quit;
