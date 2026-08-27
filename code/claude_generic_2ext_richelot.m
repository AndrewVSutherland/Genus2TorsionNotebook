// claude_generic_2ext_richelot.m — find GENERIC (End=Z) genus-2 Jacobians with
// torsion [2,22] and [2,2,14] by walking the 2-isogeny (Richelot) neighborhood of
// the abundant generic [22] and [14] curves.  A 2-isogeny preserves End^0 (so keeps
// genericity) and fixes the odd part while redistributing 2-power torsion, so a
// neighbor of a generic [22] (resp [14]) curve with 2-part (Z/2)^2 (resp (Z/2)^3)
// is a generic [2,22] (resp [2,2,14]).  Depth-3 BFS; report the first hit per target.
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);

function ModelFromEqn(fcoeffs, hcoeffs)
  f := P!fcoeffs; h := P!hcoeffs;
  return 4*f + h^2;             // y^2 = 4f+h^2 model
end function;

function IntSqfree(F)
  if F eq 0 then return F; end if;
  den := LCM([Denominator(c) : c in Coefficients(F)] cat [Integers()|1]);
  FI := [Integers()!(c*den^2) : c in Coefficients(F)];
  cont := GCD(FI);
  if cont eq 0 then return P!FI; end if;
  sq := 1; for pf in Factorization(cont) do sq *:= pf[1]^(2*(pf[2] div 2)); end for;
  return P![c div sq : c in FI];
end function;

function G2Key(F)
  return Sprintf("%o", G2Invariants(HyperellipticCurve(F)));
end function;

// End = Z certificate: root-power (geom simple) + disjoint splitting fields (center Q)
// + torsion order > 18 (not QM, LSSV). Returns bool + (p0,q0).
function EndZ(C, ordtors)
  D := Integers()!Discriminant(C); found := [];
  p0 := 0; chip := P!0;
  for p in PrimesInInterval(11, 300) do
    if D mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<a> := NumberField(chi); strict := true;
    for n in [2..12] do if Degree(MinimalPolynomial(a^n)) ne 4 then strict := false; break; end if; end for;
    dL := Degree(SplittingField(chi));
    if strict and p0 eq 0 then p0 := p; chip := chi; dP := dL; end if;
    if p0 ne 0 and p ne p0 then
      if Degree(SplittingField(chip*chi)) eq dP*dL then
        return (ordtors gt 18), p0, p;
      end if;
    end if;
    Append(~found, <p, chi, dL>);
    if #found ge 12 then break; end if;
  end for;
  return false, 0, 0;
end function;

targets := ["[2,22]", "[2,2,14]"];
seeds22 := [
  <[0,0,0,1,-2,-1,1],[1,0,1]>, <[0,1,3,4,4,3,1],[1,0,1]>, <[68,175,171,-14,-86,-21,91],[1,0,1]>,
  <[0,0,0,-4,0,4],[1,1]>, <[637,-1092,-130,429,406,-91,-143],[0,0,1]>, <[55,65,111,172,11,-35,5],[0,1,1]>,
  <[0,2,-2,1,1,-3,1],[0,1,1]>, <[0,0,0,14,-22,36,-20],[1,1]>, <[9,-9,8,-9,2,-2],[0,1,0,1]>,
  <[24,24,13,15,3,1],[0,0,1,1]>, <[63,42,-77,-224,-72,231,-63],[0,0,1]>, <[1224,-10608,14773,-5797,429,484,38],[0,0,1,1]>
];
seeds14 := [
  <[0,0,0,0,1,1],[1,0,0,1]>, <[26,37,-137,-163,296,150,-215],[1,1,1]>, <[0,0,0,0,-1],[1,0,0,1]>,
  <[-4,-4,4,4,-1,-1],[1,0,0,1]>, <[-2,24,15,24,124,-294,133],[1,1,1]>, <[0,0,1,1,2,1,1],[1,0,1]>,
  <[0,2,5,-9,-15,24,-8],[1,1]>, <[40,120,165,5,-79,32,-4],[0,0,1,1]>, <[1,-3,7,-8,8,-1],[0,1,1]>,
  <[0,1,-1,2,-1,2],[0,1,1]>, <[5,14,3,-14,-5,3,1],[1,1,1]>, <[0,0,-3,1,2,3,1],[1,0,1]>
];

procedure Hunt(seeds, want, tname, maxdepth)
  printf "==== hunting %o ====\n", tname;
  hits := 0;
  for sd in seeds do
    F0 := IntSqfree(ModelFromEqn(sd[1], sd[2]));
    if Degree(F0) lt 5 or Discriminant(F0) eq 0 then continue; end if;
    visited := {G2Key(F0)}; queue := [<F0, 0>]; cur := 1;
    while cur le #queue do
      node := queue[cur]; F := node[1]; d := node[2]; cur +:= 1;
      J := Jacobian(HyperellipticCurve(F));
      inv := Invariants(TorsionSubgroup(J));
      if inv eq want then
        ord := &*want;
        ez, p0, q0 := EndZ(HyperellipticCurve(F), ord);
        printf "CANDIDATE %o depth=%o inv=%o End=Z:%o (p0=%o,q0=%o)\n  f=%o\n", tname, d, inv, ez, p0, q0, F;
        if ez then
          printf "GENERIC_HIT %o f=%o\n", tname, F;
          hits +:= 1;
        end if;
      end if;
      if d lt maxdepth then
        Rs := [];
        try Rs := RichelotIsogenousSurfaces(J); catch e ; end try;
        for S in Rs do
          tS := Type(S);
          if tS eq JacHyp or tS eq CrvHyp then
            Cd := tS eq JacHyp select Curve(S) else S;
            fd, hd := HyperellipticPolynomials(Cd);
            FN := IntSqfree(hd eq 0 select P!fd else P!(4*fd+hd^2));
            if Degree(FN) ge 5 and Discriminant(FN) ne 0 then
              k := G2Key(FN);
              if k notin visited then Include(~visited, k); Append(~queue, <FN, d+1>); end if;
            end if;
          end if;
        end for;
      end if;
    end while;
    printf "  seed done, %o vertices\n", #queue;
    if hits ge 1 then break; end if;   // one generic example per target suffices
  end for;
  printf "%o generic hits: %o\n", tname, hits;
end procedure;

Hunt(seeds22, [2,22], "[2,22]", 3);
Hunt(seeds14, [2,2,14], "[2,2,14]", 3);
print "ALL_DONE";
quit;
