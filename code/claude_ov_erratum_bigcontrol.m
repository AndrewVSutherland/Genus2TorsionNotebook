// claude_ov_erratum_bigcontrol.m -- genus2red ERRATUM propagation (2026-07-25).
//
// THE BIG CONTROL, aimed squarely at the one soft spot in the replacement
// method.  Magma's Conductor for a genus-2 curve prints
//     "WARNING: Using Ogg's formula when v_2(D)>=12, no correctness guarantee"
// whenever v_2(disc_min) >= 12, and four of the eleven contact-7 [2,2,14]
// curves (plus the BLP C4 [2,22] curve) sit in exactly that regime.  So: how
// good is the Ogg fallback there, empirically?
//
// Control set: EVERY curve in the production LMFDB g2c_curves table with
// 2^12 | abs_disc (1113 of them, pulled 2026-07-25 into
// data/claude_ov_erratum_g2c_v2disc12.csv, converted to
// data/claude_ov_erratum_g2c_control.m).  Their conductors are published data
// (Booker-Sijsling-Sutherland-Voight), independent of this project.  v_2(cond)
// runs from 0 to 19 over the set.
//
// Forked, 12 children (raw POSIX fork; see code/claude_fork_template.m).
//
// Usage:   magma -b NCH:=12 MemGB:=4 code/claude_ov_erratum_bigcontrol.m
// Markers: BIGCTRL (mismatches + Ogg-regime rows) / BIGSUMMARY / ERRATUM_BIGCONTROL_DONE
SetColumns(0);
if not assigned NCH   then NCH := 12;   elif Type(NCH)   eq MonStgElt then NCH   := StringToInteger(NCH);   end if;
if not assigned MemGB then MemGB := 4;  elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

load "data/claude_ov_erratum_g2c_control.m";     // defines BIGCTRL
n := #BIGCTRL;
printf "BIGCTRL loaded %o curves, forking %o children\n", n, NCH;

dir := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/erratum_bigctrl";
System("mkdir -p " cat dir);

for c in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then                                        // ---- child ----
    F := Open(dir cat "/part" cat IntegerToString(c) cat ".txt", "w");
    for i in [c+1..n by NCH] do
      r := BIGCTRL[i];
      lab := r[1]; Npub := r[2]; adisc := r[3];
      fc := r[4]; hc := r[5];
      ok := true; msg := "";
      try
        f := P!fc; h := P!hc;
        C := HyperellipticCurve(f, h);
        Cm := ReducedMinimalWeierstrassModel(C);
        Dm := Z!Discriminant(Cm);
        v2D := Valuation(AbsoluteValue(Dm), 2);
        N := Conductor(Cm);
        Puts(F, lab cat " " cat IntegerToString(Npub) cat " " cat IntegerToString(N)
                cat " " cat IntegerToString(v2D) cat " " cat IntegerToString(adisc));
      catch ee
        Puts(F, lab cat " " cat IntegerToString(Npub) cat " ERR -1 " cat IntegerToString(adisc));
      end try;
    end for;
    Flush(F); delete F;
    quit;                                                 // child MUST exit
  end if;
end for;
WaitForAllChildren();                                     // ---- parent ----

nall := 0; nerr := 0; nmatch := 0; nmis := 0;
noggall := 0; noggmatch := 0;
mis := []; v2set := {};
for c in [0..NCH-1] do
  for l in Split(Read(dir cat "/part" cat IntegerToString(c) cat ".txt"), "\n") do
    if #l eq 0 then continue; end if;
    s := Split(l, " ");
    lab := s[1]; Npub := StringToInteger(s[2]);
    if s[3] eq "ERR" then
      nerr +:= 1; printf "BIGCTRL ERROR %o\n", lab; continue;
    end if;
    N := StringToInteger(s[3]); v2D := StringToInteger(s[4]);
    nall +:= 1;
    Include(~v2set, Valuation(Npub,2));
    if v2D ge 12 then
      noggall +:= 1;
      if N eq Npub then noggmatch +:= 1; end if;
    end if;
    if N eq Npub then nmatch +:= 1; else
      nmis +:= 1; Append(~mis, <lab, Npub, N, v2D>);
      printf "BIGCTRL MISMATCH %o  Npub=%o  Nmagma=%o  v_2(disc_min)=%o\n", lab, Npub, N, v2D;
    end if;
  end for;
end for;
printf "BIGSUMMARY curves=%o errors=%o matched=%o mismatched=%o\n", nall, nerr, nmatch, nmis;
printf "BIGSUMMARY Ogg-warned regime (v_2(disc_min)>=12): %o curves, %o matched, %o mismatched\n",
    noggall, noggmatch, noggall - noggmatch;
printf "BIGSUMMARY v_2(published conductor) values seen: %o\n", Sort(Setseq(v2set));
print "ERRATUM_BIGCONTROL_DONE";
quit;
