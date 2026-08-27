SetColumns(0); Q:=Rationals(); P<x>:=PolynomialRing(Q);
// How often does a Richelot isogeny from a 2-rank-1 genus-2 Jacobian RAISE the 2-rank?
// Sample: f = q * Q4 with q a rational quadratic and Q4 an irreducible quartic whose
// resolvent cubic has a rational root (=> Galois-stable, non-split, (2,2)-kernel).
cnt := 0; raised := 0; tot := 0;
for a in [-4..4] do for b in [-4..4] do for c in [-3..3] do
  if cnt ge 400 then break; end if;
  Q4 := x^4 + a*x^3 + b*x^2 + c*x + 1;
  if not IsIrreducible(Q4) then continue; end if;
  G := GaloisGroup(Q4);
  if #G gt 8 then continue; end if;                 // Galois subset of D4
  q := x^2 + x + 1;
  f := q*Q4;
  if Degree(f) ne 6 or Discriminant(f) eq 0 then continue; end if;
  C := HyperellipticCurve(f);
  if not IsIrreducible(f) then end if;
  J := Jacobian(C);
  r0 := #Invariants(TwoTorsionSubgroup(J));
  if r0 ne 1 then continue; end if;                 // want 2-rank exactly 1
  cnt +:= 1;
  RS := RichelotIsogenousSurfaces(J);
  nb := 0;
  for S in RS do
    if Type(S) ne JacHyp then continue; end if;
    nb +:= 1; tot +:= 1;
    r1 := #Invariants(TwoTorsionSubgroup(S));
    if r1 gt r0 then raised +:= 1; end if;
  end for;
  printf "Q4=%o Gal=%o  #neighbours=%o  ranks=%o\n", Q4, #G, nb,
     [#Invariants(TwoTorsionSubgroup(S)) : S in RS | Type(S) eq JacHyp];
end for; end for; end for;
printf "SAMPLED %o rank-1 curves, %o neighbours, %o with HIGHER 2-rank\n", cnt, tot, raised;
quit;
