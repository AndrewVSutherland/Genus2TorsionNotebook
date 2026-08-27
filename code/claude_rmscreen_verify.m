// claude_rmscreen_verify.m -- independent reproduction of the broken RM screen
// and its one-line repair (orchestrator check of Lane 3's finding, 2026-07-25).
//
// The lab doctrine "real-subfield-disc census CONSTANT = RM, SCATTERING = End=Z"
// is false in the second direction. At a prime inert in the field of definition
// of the RM, a1 = Tr(a_p) = 0 exactly, so n_p = a1^2-4(a2-2p) = 8p-4a2 and its
// squarefree core is junk varying with p. Such primes have density 1/2, so a
// geometrically-RM curve whose RM is NOT defined over Q scatters and is
// misread as End = Z.
//
// REPAIR: count only primes with a1 <> 0 and n_p neither 0 nor a perfect square.
//
// Test curve: LMFDB 12500.b.50000.1, certified geom_end_alg = 'RM' with
// is_gl2_type = FALSE (i.e. exactly the failing class).
// Result (results/claude_rmscreen_verify.log):
//     primes used 60, a1 = 0 at 33 of them (55%)
//     NAIVE    -> 23 distinct cores -> verdict "End=Z"   WRONG
//     REPAIRED ->  1 distinct core {5} -> verdict "RM"   CORRECT (and names sqrt5)
//
// Does NOT affect: the constant-census => RM direction, nor any strict two-prime
// Frobenius certificate (a theorem; an even chi fails root-power strictness and
// is refused rather than wrongly certified).

SetColumns(0); SetMemoryLimit(4*10^9);
Z:=Integers(); Q:=Rationals(); P<x>:=PolynomialRing(Q);
// LMFDB 12500.b.50000.1 : geom_end_alg = RM, is_gl2_type = FALSE
f := P![-2,3,0,-3,0,1];  h := P![1,0,0,1];
C := HyperellipticCurve(f,h);
F := 4*f + h^2; Cs := HyperellipticCurve(F);
D := Z!Discriminant(Cs);
naive := {}; repaired := {}; na1z := 0; np := 0;
for p in PrimesInInterval(3,300) do
  if D mod p eq 0 then continue; end if;
  chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(Cs,GF(p))))));
  a1 := -Coefficient(chi,3); a2 := Coefficient(chi,2);
  n := a1^2 - 4*(a2 - 2*p);
  np +:= 1;
  if a1 eq 0 then na1z +:= 1; end if;
  if n ne 0 then Include(~naive, Squarefree(Z!Abs(n))); end if;
  // REPAIR: only primes with a1 <> 0 and n neither 0 nor a perfect square
  if a1 ne 0 and n ne 0 and not IsSquare(Z!Abs(n)) then Include(~repaired, Squarefree(Z!Abs(n))); end if;
end for;
printf "primes used = %o ; primes with a1 = 0 : %o  (%o%%)\n", np, na1z, Round(100*na1z/np);
printf "NAIVE screen    -> %o distinct cores %o  => verdict %o\n", #naive,
   Sort(Setseq(naive))[1..Min(8,#naive)], (#naive le 1) select "RM" else "End=Z (WRONG)";
printf "REPAIRED screen -> %o distinct cores %o  => verdict %o\n", #repaired,
   Sort(Setseq(repaired)), (#repaired le 1) select "RM (CORRECT)" else "End=Z";
quit;
