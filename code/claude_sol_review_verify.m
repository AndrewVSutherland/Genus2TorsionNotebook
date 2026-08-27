// claude_sol_review_verify.m — verify the load-bearing numeric claims in ChatGPT_summary.tex
// and settle the certificate strength of the cyclic [40] and [28] curves
// (notes record only Q-simplicity; the summary claims geometric simplicity).
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

// strict root-power test at p: chi irreducible AND [Q(pi^n):Q]=4 for n=2..12
function StrictAt(C, p)
  Cp := ChangeRing(C, GF(p));
  chi := Px!Reverse(Coefficients(EulerFactor(Jacobian(Cp))));
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false, chi; end if;
  K<pi> := NumberField(chi);
  for n in [2..12] do
    if Degree(MinimalPolynomial(pi^n)) ne 4 then return false, chi; end if;
  end for;
  return true, chi;
end function;

procedure StrictSweep(C, label, plo, phi, want)
  D := Integers()!Discriminant(HyperellipticPolynomials(C));
  good := 0;
  for p in PrimesInInterval(plo, phi) do
    if good ge want then break; end if;
    if D mod p eq 0 then continue; end if;
    ok, chi := StrictAt(C, p);
    if ok then printf "%o: STRICT pass at p=%o\n", label, p; good +:= 1; end if;
  end for;
  printf "%o: strict passes found: %o\n", label, good;
end procedure;

print "==== cyclic [40] curve (contact5 t=-1/3) ====";
f40 := -324*x^5 + 1296*x^4 + 1944*x^3 - 5103*x^2 - 4374*x + 6561;
C40 := HyperellipticCurve(f40);
T40 := TorsionSubgroup(Jacobian(C40));
printf "[40] torsion: %o (order %o)\n", Invariants(T40), #T40;
StrictSweep(C40, "[40]", 7, 300, 2);

print "==== cyclic [28] curve (contact7 rational root) ====";
f28 := 4*x^5 + 21*x^4 - 70*x^3 + 79*x^2 - 42*x + 9;
C28 := HyperellipticCurve(f28);
T28 := TorsionSubgroup(Jacobian(C28));
printf "[28] torsion: %o (order %o)\n", Invariants(T28), #T28;
// confirm p=5 is structurally useless (a_p = 0)
_, chi5 := StrictAt(C28, 5);
printf "[28] chi at p=5: %o\n", chi5;
StrictSweep(C28, "[28]", 3, 300, 2);

print "==== record curve: Sol's reduction and Frobenius claims ====";
f1 := 756900*x^6 + 737595570*x^5 + 150572203590*x^4 - 15854483576121*x^3 - 530648977741620*x^2 + 32014154874551031*x + 830742747091037849;
Crec := HyperellipticCurve(f1, x^2+1);
Csim := SimplifiedModel(Crec);
n31 := #Jacobian(ChangeRing(Csim, GF(31)));
n37 := #Jacobian(ChangeRing(Csim, GF(37)));
printf "#J(F_31) = %o (claim 864), #J(F_37) = %o (claim 1248), gcd = %o\n", n31, n37, GCD(n31, n37);
for p in [37, 73] do
  chi := Px!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(Csim, GF(p))))));
  printf "P_%o = %o\n", p, chi;
end for;

print "==== C_66 claims ====";
f66 := 1872*x^5 - 3000*x^4 + 6969*x^3 - 1691*x^2 + 4875*x;
C66 := HyperellipticCurve(f66);
T66 := TorsionSubgroup(Jacobian(C66));
printf "C66 torsion: %o\n", Invariants(T66);
printf "#J(F_7) = %o (claim 36), #J(F_11) = %o (claim 144)\n",
  #Jacobian(ChangeRing(C66, GF(7))), #Jacobian(ChangeRing(C66, GF(11)));

print "ALL_DONE";
quit;
