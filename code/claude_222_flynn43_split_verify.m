// claude_222_flynn43_split_verify.m — VERIFY the [2,22] candidates from the
// quadratic-splitting cover of the parametrized Flynn (4,3) QF family.
// The cover Y^2 = 117 tau^4 - 3673 tau^3 + 40642 tau^2 - 168108 tau + 107272 is
// genus 1 with rational points tau = -14, 69/4, 8446/711, 3981/5171 (and tau=11
// on the boundary Y=0).  At such tau, disc(x^2+u x+v) is a square, so the
// Flynn member F_{t(tau)} has factor type [1,1,1,3] => 2-rank 2; with the
// order-11 infinity class: torsion should contain [2,22].
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);
FF<T> := FunctionField(Q);

ufun := (-11/100*T^3 + 919/50*T^2 - 8017/25*T + 40154/25)/(T^3 - 23*T^2 + 168*T - 396);
vfun := (2401/1600*T^4 - 12593/200*T^3 + 205647/200*T^2 - 385757/50*T + 2253001/100)/(T^4 - 34*T^3 + 421*T^2 - 2244*T + 4356);
tfun := (-21609/40000*T^6 + 435561/10000*T^5 - 2897123/2000*T^4 + 6390763/250*T^3 - 126932427/500*T^2 + 844854361/625*T - 1894773841/625)/(T^6 - 46*T^5 + 865*T^4 - 8520*T^3 + 46440*T^2 - 133056*T + 156816);

function EndZ(C, ord)
  D := Integers()!Discriminant(C); p0 := 0; chip := P!0; dP := 0; nf := 0;
  for p in PrimesInInterval(11, 400) do
    if D mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<aa> := NumberField(chi); st := true;
    for n in [2..12] do if Degree(MinimalPolynomial(aa^n)) ne 4 then st := false; break; end if; end for;
    dL := Degree(SplittingField(chi));
    if st and p0 eq 0 then p0 := p; chip := chi; dP := dL; continue; end if;
    if p0 ne 0 and st and Degree(SplittingField(chip*chi)) eq dP*dL then return (ord gt 18), p0, p, chip, chi; end if;
    nf +:= 1; if nf ge 20 then break; end if;
  end for;
  return false, 0, 0, P!0, P!0;
end function;

taus := [ -14, 69/4, 8446/711, 3981/5171 ];
for tau in taus do
  u0 := Evaluate(ufun, tau); v0 := Evaluate(vfun, tau); t0 := Evaluate(tfun, tau);
  F := x^6+2*x^5+(2*t0+3)*x^4+2*x^3+(t0^2+1)*x^2+2*t0*(1-t0)*x+t0^2;
  if Discriminant(F) eq 0 then printf "tau=%o DEGENERATE\n", tau; continue; end if;
  // integral model
  den := LCM([Denominator(co) : co in Coefficients(F)]);
  Fi := P![ co*den^2 : co in Coefficients(F) ];
  ftype := Sort([Degree(pe[1]) : pe in Factorization(Fi)]);
  C := HyperellipticCurve(Fi);
  J := Jacobian(C);
  inv := Invariants(TorsionSubgroup(J));
  printf "TAU=%o t=%o type=%o TORSION=%o\n", tau, t0, ftype, inv;
  if inv eq [2,22] then
    ez, p0, q0, cp, cq := EndZ(C, 44);
    printf "  EndZ=%o p0=%o q0=%o\n", ez, p0, q0;
    if ez then
      printf "GENERIC_222_VERIFIED tau=%o\n  chi%o=%o\n  chi%o=%o\n  F=%o\n", tau, p0, cp, q0, cq, F;
    end if;
  end if;
end for;
print "SPLIT_VERIFY_DONE";
quit;
