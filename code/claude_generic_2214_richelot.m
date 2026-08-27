// claude_generic_2214_richelot.m — generic [2,2,14] by Richelot from generic [2,14].
// A generic [2,14] curve has 2-rank 2, so admits rational (2,2)-isogenies; a Richelot
// codomain has the same End^0 (generic) and 11->7 odd part fixed, and can carry 2-rank 3
// => generic [2,2,14].  Depth-4 BFS over the 2-isogeny neighborhood of each seed.
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);

function IntSqfree(F)
  if F eq 0 then return F; end if;
  den := LCM([Denominator(c) : c in Coefficients(F)] cat [Integers()|1]);
  FI := [Integers()!(c*den^2) : c in Coefficients(F)];
  cont := GCD(FI); if cont eq 0 then return P!FI; end if;
  sq := 1; for pf in Factorization(cont) do sq *:= pf[1]^(2*(pf[2] div 2)); end for;
  return P![c div sq : c in FI];
end function;
function ModelFromEqn(fc, hc) return 4*(P!fc) + (P!hc)^2; end function;
function G2Key(F) return Sprintf("%o", G2Invariants(HyperellipticCurve(F))); end function;

function EndZ(C, ord)
  D := Integers()!Discriminant(C); p0 := 0; chip := P!0; dP := 0; nfound := 0;
  for p in PrimesInInterval(11, 300) do
    if D mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<a> := NumberField(chi); strict := true;
    for n in [2..12] do if Degree(MinimalPolynomial(a^n)) ne 4 then strict := false; break; end if; end for;
    dL := Degree(SplittingField(chi));
    if strict and p0 eq 0 then p0 := p; chip := chi; dP := dL; continue; end if;
    if p0 ne 0 then
      if Degree(SplittingField(chip*chi)) eq dP*dL then return (ord gt 18), p0, p, chip, chi; end if;
    end if;
    nfound +:= 1; if nfound ge 14 then break; end if;
  end for;
  return false, 0, 0, P!0, P!0;
end function;

seeds := [
  <[0,0,-2,1,0,-1,1],[1,0,1]>, <[0,-3,9,13,-3,-3,1],[0,1,1]>, <[9,-3,0,-3,-4],[0,1,0,1]>,
  <[0,-6,-21,26,66,-102,36],[0,1,1]>, <[0,12,-45,21,36,11,1],[0,0,1]>, <[1,-1,0,-3,-8,0,2],[0,1,0,1]>,
  <[-34,-135,172,360,486,35,-5],[1,0,1]>, <[1,-11,30,-3,-20,0,2],[0,1,0,1]>, <[96,-249,68,707,91,168,63],[1,1]>,
  <[0,-12,13,-15,-3,1],[0,0,1,1]>, <[0,12,9,4,3,-3,1],[0,1,1]>, <[0,1,0,-6,-2,-6],[1,1]>,
  <[0,1,6,6,34,-48],[1,1]>, <[-35,465,-1554,510,2395,960,105],[0,1]>
];

maxdepth := 4; hits := 0;
for si in [1..#seeds] do
  sd := seeds[si];
  F0 := IntSqfree(ModelFromEqn(sd[1], sd[2]));
  if Degree(F0) lt 5 or Discriminant(F0) eq 0 then continue; end if;
  visited := {G2Key(F0)}; queue := [<F0,0>]; cur := 1; nrich := 0;
  while cur le #queue do
    node := queue[cur]; F := node[1]; d := node[2]; cur +:= 1;
    J := Jacobian(HyperellipticCurve(F));
    inv := Invariants(TorsionSubgroup(J));
    if inv eq [2,2,14] then
      ez, p0, q0, cp, cq := EndZ(HyperellipticCurve(F), 56);
      printf "CANDIDATE seed%o depth=%o inv=%o End=Z:%o (p0=%o,q0=%o)\n  f=%o\n", si, d, inv, ez, p0, q0, F;
      if ez then
        printf "GENERIC_HIT_2214 f=%o\n  chi_p0(%o)=%o  chi_q0(%o)=%o\n", F, p0, cp, q0, cq;
        hits +:= 1; break;
      end if;
    end if;
    if d lt maxdepth then
      Rs := [];
      try Rs := RichelotIsogenousSurfaces(J); catch e ; end try;
      for S in Rs do
        tS := Type(S);
        if tS eq JacHyp or tS eq CrvHyp then
          nrich +:= 1;
          Cd := tS eq JacHyp select Curve(S) else S;
          fd, hd := HyperellipticPolynomials(Cd);
          FN := IntSqfree(hd eq 0 select P!fd else P!(4*fd+hd^2));
          if Degree(FN) ge 5 and Discriminant(FN) ne 0 then
            k := G2Key(FN);
            if k notin visited then Include(~visited, k); Append(~queue, <FN,d+1>); end if;
          end if;
        end if;
      end for;
    end if;
  end while;
  printf "seed%o: %o vertices, %o richelot-codomains\n", si, #queue, nrich;
  if hits ge 1 then break; end if;
end for;
printf "[2,2,14] generic hits: %o\n", hits;
print "ALL_DONE";
quit;
