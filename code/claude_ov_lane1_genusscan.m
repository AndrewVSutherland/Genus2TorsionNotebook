// claude_ov_lane1_genusscan.m -- Lane 1: genus of the u = c slice C_c of the
// contact-7 three-root surface R(s,t,u)=0, for every rational c of bounded
// height.  Generic C_c is a (3,3) curve of arithmetic genus 4; c = -1/2 is
// known to drop to genus 1 (rank 1).  Which other c drop?
// Shards by Part/NParts.  Usage: magma -b H:=40 NParts:=192 Part:=k ...
SetColumns(0);
if not assigned H then H := 20; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned MemGB then MemGB := 2; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Rst<s,t> := PolynomialRing(Q,2);
A2 := AffineSpace(Rst);

// R(s,t,u) with u specialised to c, written out explicitly (verified against the
// gp rebuild in results/claude_ov_lane1_surface.log).
function Rslice(c)
  c3 := (2*t^2+4*t+2)*s^3 + (2*t^3+8*t^2+10*t+4)*s^2 + (4*t^3+10*t^2+8*t+2)*s + (2*t^3+4*t^2+2*t);
  c2 := (2*t^3+8*t^2+10*t+4)*s^3 + (8*t^3+20*t^2+16*t+4)*s^2 + (10*t^3+16*t^2+8*t+1)*s + (4*t^3+4*t^2+t);
  c1 := (4*t^3+10*t^2+8*t+2)*s^3 + (10*t^3+16*t^2+8*t+1)*s^2 + (8*t^3+8*t^2+2*t)*s + (2*t^3+t^2);
  c0 := (2*t^3+4*t^2+2*t)*s^3 + (4*t^3+4*t^2+t)*s^2 + (2*t^3+t^2)*s;
  return c3*c^3 + c2*c^2 + c1*c + c0;
end function;

vals := [];
for d in [1..H] do for nu in [-H..H] do
  if nu ne 0 and Gcd(Abs(nu),d) eq 1 then Append(~vals, Q!nu/d); end if;
end for; end for;
vals := Sort(Setseq(Seqset(vals)));
if Part eq 0 then printf "GENUSSCAN H=%o #c-values=%o NParts=%o\n", H, #vals, NParts; end if;

cnt := 0;
for i in [1..#vals] do
  if (i-1) mod NParts ne Part then continue; end if;
  c := vals[i];
  if c eq 0 or c eq -1 then continue; end if;   // degenerate chart values
  F := Rslice(c);
  if F eq 0 then printf "IDENTICALLY_ZERO c=%o\n", c; continue; end if;
  ok := true; g := -1;
  try
    Cv := Curve(A2, F);
    if IsIrreducible(Cv) then g := Genus(Cv); else g := -2; end if;
  catch e ok := false; end try;
  cnt +:= 1;
  if not ok then printf "ERROR c=%o\n", c; continue; end if;
  if g ne 4 then printf "LOWGENUS c=%o genus=%o\n", c, g; end if;
  if cnt mod 200 eq 0 then printf "PROGRESS part=%o done=%o at c=%o\n", Part, cnt, c; end if;
end for;
printf "GENUSSCAN_DONE part=%o count=%o\n", Part, cnt;
quit;
