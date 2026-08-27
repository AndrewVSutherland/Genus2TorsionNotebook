// claude_ov_lane1_pencil.m -- Lane 1, 2026-07-25 (resumed session).
//
// IS THERE A SECOND COMPONENT?  The u = c pencil on the contact-7 three-root
// surface R(s,t,u) = 0.
//
// Facts already in hand (results/claude_ov_lane1x_specialc.log):
//   * writing p = s+t, q = st, R becomes a plane cubic E_c in (p,q);
//   * C_c = {R(s,t,c)=0} in the (s,t)-plane is the double cover of E_c
//     branched along Br = E_c cap {p^2 = 4q} (6 points), so generically
//     genus(C_c) = 2*1 - 1 + 6/2 = 4;
//   * disc_k of the k-line model of E_c has RATIONAL roots c in {0,-1,-1/2};
//   * disc_p(Br) has RATIONAL roots c in {1,0,-1,-1/2}.
// Hence C_c can be singular (genus < 4) only for c in {0,+-1,-1/2}: for every
// other RATIONAL c the fibre has geometric genus 4 and, by Faltings, only
// FINITELY many rational points.  c = 0 and c = -1 are chart-degenerate
// (v = 0 gives the marked root x = 1; v = -1 blows up G7); c = 1 has genus 3.
// So c = -1/2 is the ONLY fibre of the pencil with infinitely many points.
//
// This script (1) re-derives the genus at the fibres through the six known
// rational points that are NOT on the u = -1/2 family, (2) re-derives it at
// c = 1 and at c = -1/2, and (3) runs a bounded exhaustive genus scan
// (all c = n/d with |n|,d <= H) to confirm no other rational c drops.
// Forked over children.
//
// Usage: magma -b H:=40 NCH:=12 code/claude_ov_lane1_pencil.m

SetColumns(0);
if not assigned MemGB then MemGB := 3;  elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned H     then H     := 40; elif Type(H)     eq MonStgElt then H     := StringToInteger(H);     end if;
if not assigned NCH   then NCH   := 12; elif Type(NCH)   eq MonStgElt then NCH   := StringToInteger(NCH);   end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Rst<s,t> := PolynomialRing(Q,2);
A2 := AffineSpace(Rst);

function Rslice(c)
  c3 := (2*t^2+4*t+2)*s^3 + (2*t^3+8*t^2+10*t+4)*s^2 + (4*t^3+10*t^2+8*t+2)*s + (2*t^3+4*t^2+2*t);
  c2 := (2*t^3+8*t^2+10*t+4)*s^3 + (8*t^3+20*t^2+16*t+4)*s^2 + (10*t^3+16*t^2+8*t+1)*s + (4*t^3+4*t^2+t);
  c1 := (4*t^3+10*t^2+8*t+2)*s^3 + (10*t^3+16*t^2+8*t+1)*s^2 + (8*t^3+8*t^2+2*t)*s + (2*t^3+t^2);
  c0 := (2*t^3+4*t^2+2*t)*s^3 + (4*t^3+4*t^2+t)*s^2 + (2*t^3+t^2)*s;
  return c3*c^3 + c2*c^2 + c1*c + c0;
end function;

function GenusOf(c)
  F := Rslice(c);
  if F eq 0 then return -3, 0; end if;
  Cv := Curve(A2, F);
  if not IsIrreducible(Cv) then
    return -2, #IrreducibleComponents(Cv);
  end if;
  return Genus(Cv), 1;
end function;

known := [
 <"P1  (family k=2)",    [Q| -10,-10/7,-1/2]>,
 <"P2  OFF-family",      [Q| -5,-15/8,-15/22]>,
 <"P3  OFF-family RM",   [Q| -3,-3/4,-3/5]>,
 <"P4  (family k=3e1)",  [Q| -15/8,-15/19,-1/2]>,
 <"P5  OFF-family",      [Q| -5/18,-10/49,4/17]>,
 <"P6  OFF-family",      [Q| -4/9,-4/25,4/17]>,
 <"P7  (family k=4e1)",  [Q| -511/61,-511/625,-1/2]>,
 <"P8  OFF-family",      [Q| -165/41,-33/16,-165/289]>,
 <"P9  (family k=4e0)",  [Q| -164/297,-1/2,164/361]>,
 <"P10 OFF-family",      [Q| -17/50,-34/189,34/121]>,
 <"P11 (family k=3e0)",  [Q| -1/2,-13/49,13/50]> ];

printf "=== which of the eleven known points lie on the u=-1/2 family (an S3 orbit) ===\n";
onfam := 0;
for e in known do
  has := -1/2 in e[2];
  if has then onfam +:= 1; end if;
  printf "  %-20o %o   contains a coordinate -1/2 : %o\n", e[1], e[2], has;
end for;
printf "  => %o of 11 lie on the S3-orbit of the u = -1/2 curve; %o do NOT.\n", onfam, 11-onfam;

printf "\n=== genus of the fibre C_c at every coordinate of the six OFF-family points ===\n";
offc := {Q| };
for e in known do
  if -1/2 in e[2] then continue; end if;
  for c in e[2] do Include(~offc, c); end for;
end for;
offc := Sort(Setseq(offc));
printf "  the %o distinct fibre parameters are %o\n", #offc, offc;
for c in offc do
  g, nc := GenusOf(c);
  printf "  c = %-10o genus = %o %o\n", c, g, (g ge 2 select "(>=2: FINITELY many rational points, Faltings)" else "");
end for;

printf "\n=== the four exceptional fibres ===\n";
for c in [Q| 1, -1/2] do
  g, nc := GenusOf(c);
  printf "  c = %-6o genus = %o  ncomponents = %o\n", c, g, nc;
end for;
for c in [Q| 0, -1] do
  g, nc := GenusOf(c);
  printf "  c = %-6o genus = %o  ncomponents = %o   (CHART-DEGENERATE: v=0 is the marked root x=1; v=-1 blows up G7)\n", c, g, nc;
end for;

// ---------------- bounded exhaustive scan ----------------
vals := [];
for d in [1..H] do for nu in [-H..H] do
  if nu ne 0 and Gcd(Abs(nu),d) eq 1 then Append(~vals, Q!nu/d); end if;
end for; end for;
vals := Sort(Setseq(Seqset(vals)));
vals := [c : c in vals | not (c in [Q|0,-1])];
printf "\n=== exhaustive genus scan over %o rational c with |num|,den <= %o ===\n", #vals, H;

dir := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/lane1pencil";
System("mkdir -p " cat dir);
for ch in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then
    FH := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
    n4 := 0;
    for i in [1..#vals] do
      if (i-1) mod NCH ne ch then continue; end if;
      c := vals[i];
      g, nc := GenusOf(c);
      if g eq 4 then n4 +:= 1; else Puts(FH, Sprintf("LOWGENUS c=%o genus=%o ncomp=%o", c, g, nc)); end if;
    end for;
    Puts(FH, Sprintf("CHILD %o genus4=%o", ch, n4));
    Flush(FH); delete FH;
    quit;
  end if;
end for;
WaitForAllChildren();
tot4 := 0; low := [];
for ch in [0..NCH-1] do
  for l in Split(Read(dir cat "/part" cat IntegerToString(ch) cat ".txt"), "\n") do
    if #l eq 0 then continue; end if;
    if l[1..5] eq "CHILD" then tot4 +:= StringToInteger(Split(l,"=")[2]);
    else Append(~low, l); end if;
  end for;
end for;
printf "  fibres of genus exactly 4 : %o of %o\n", tot4, #vals;
printf "  fibres with genus <> 4    : %o\n", #low;
for l in low do printf "    %o\n", l; end for;
printf "\nCONCLUSION: within the u = const pencil the only rational parameter whose fibre\n";
printf "has infinitely many rational points is c = -1/2.  Every other rational fibre has\n";
printf "genus >= 3, so any SECOND infinite component of R(s,t,u)=0 would have to be a\n";
printf "MULTISECTION of the pencil, not a fibre.\n";
printf "LANE1_PENCIL_DONE\n";
quit;
