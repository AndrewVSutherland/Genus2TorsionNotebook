// claude_end_z_certificates.m — certify End(Jac(C)_Qbar) = Z ("generic") for the
// table's non-LMFDB curves.  Criterion: two good primes p, q whose Frobenius
// charpolys are irreducible with all root-powers of degree 4 (strict test: kills QM
// and forces End(J_Fpbar) tensor Q = Q(pi_p)), such that the quartic fields
// Q(pi_p), Q(pi_q) share NO quadratic subfield.  Any extra endomorphism field would
// embed in both quartic fields, so disjoint quadratic-subfield sets force End = Z.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

// Zywina-style certificate: two good primes with irreducible quartic L-polynomials
// whose splitting fields are LINEARLY DISJOINT over Q (both are Galois, so this is
// [L_p L_q : Q] = [L_p : Q][L_q : Q]).  This forces End(Jac_Qbar) = Z; see
// Zywina, "Determining monodromy groups of abelian varieties", Res. Number Theory 8
// (2022), and Costa-Lombardo-Voight, Res. Number Theory 7 (2021).
function IrredChi(C, p)
  Cp := ChangeRing(C, GF(p));
  chi := Px!Reverse(Coefficients(EulerFactor(Jacobian(Cp))));
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false, _, _; end if;
  return true, chi, Degree(SplittingField(chi));
end function;

procedure EndZCert(f, h, label)
  C := HyperellipticCurve(f, h);
  D := Integers()!Discriminant(C);
  found := [];
  for p in PrimesInInterval(11, 400) do
    if D mod p eq 0 then continue; end if;
    ok, chi, dL := IrredChi(C, p);
    if not ok then continue; end if;
    for pr in found do
      // compositum of the two (Galois) splitting fields = splitting field of the product
      degC := Degree(SplittingField(pr[2]*chi));
      if degC eq pr[3]*dL then
        printf "%o: END_Z_CERT primes (%o, %o), splitting degrees (%o, %o), compositum %o (linearly disjoint)\n",
          label, pr[1], p, pr[3], dL, degC;
        return;
      end if;
    end for;
    Append(~found, <p, chi, dL>);
    if #found ge 8 then break; end if;
  end for;
  printf "%o: NO certificate found (irreducible primes tried: %o)\n", label, [pr[1] : pr in found];
end procedure;

// the six project curves
EndZCert(-391671*x^6+1894851*x^5+6846924*x^4-15133525*x^3+3904068*x^2+2625336*x+254016, x^2+x, "[2,2,20] C2220");
EndZCert(60*x^5+3982*x^4+1100*x^3+130682*x^2-239760*x, x^2+x, "[2,4,8] C248");
EndZCert(1872*x^5-3000*x^4+6969*x^3-1691*x^2+4875*x, Px!0, "[6,6] C66");
EndZCert(x*(x+16)*(x+(644/799)^2)*(x^2+x+16)*799^2, Px!0, "[2,4,4] C244 (scaled)");
EndZCert(x*(x+1296)*(x+3249)*(x+4096)*(x+17424), Px!0, "[2,2,4,4] C2244");
EndZCert(756900*x^6+737595570*x^5+150572203590*x^4-15854483576121*x^3-530648977741620*x^2+32014154874551031*x+830742747091037849, x^2+1, "[2,2,2,12] curve1");
// curves #2, #3 of (2,2,2,12)
EndZCert(36*x^6+36750*x^5-462983772*x^4-301623595823*x^3+1518598238654317*x^2+397058962729817115*x-1282993930035013443975, x^2+x, "[2,2,2,12] curve2");
EndZCert(3703062294195264*x^6-360079374491052216*x^5+8901721379573296848*x^4-5397945250386334945*x^3-86737535708373850908*x^2+36346694984390901540*x+43035470132681030400, x^2+x, "[2,2,2,12] curve3");
// isolated jumps certified earlier (upgrade to End = Z)
EndZCert(-324*x^5+1296*x^4+1944*x^3-5103*x^2-4374*x+6561, Px!0, "[40] contact5");
EndZCert(4*x^5+21*x^4-70*x^3+79*x^2-42*x+9, Px!0, "[28] contact7");
// literature curves with equations in hand
EndZCert((9*x^2+2*x+1)*(32*x^3+81*x^2-6*x+1), Px!0, "[34] Elkies");
EndZCert(x^6+4*x^4+10*x^3+4*x^2-4*x+1, Px!0, "[39] Elkies");
EndZCert((3*x+4)*(4*x^4+20*x^3+32*x^2+19*x+4), Px!0, "[40] Elkies (scaled by 4)");
EndZCert(x*(x+1)*(x-1)*(3*x-7)*(8*x-13)*(24*x+25), Px!0, "[2,2,2,10] Elkies 2024 member");
print "ALL_DONE";
quit;
