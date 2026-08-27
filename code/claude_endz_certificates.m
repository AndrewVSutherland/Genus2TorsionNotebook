// claude_endz_certificates.m — publication certificates for the displayed
// non-database curves of paper/torsion_realizations.tex.  For each curve emits, to
// data/claude_endz_certificates.txt, a self-contained End(Jac_Qbar) = Z certificate:
//   (i)   exact TorsionSubgroup (the displayed invariant factors);
//   (ii)  a good prime p_0 with L-polynomial chi_{p_0} irreducible AND every root
//         power pi^n (2<=n<=12) of degree 4  =>  Jac is geometrically simple
//         (the Leprevost/root-power criterion; simplicity-certificates skill);
//   (iii) a second good prime q_0 whose L-polynomial chi_{q_0} is irreducible with
//         splitting field linearly disjoint from that of chi_{p_0}  =>  the center of
//         End(Jac_Qbar) is Q, excluding real-multiplication and CM (Costa-Lombardo-
//         Voight 2021; Zywina 2022);
//   (iv)  #Jac(Q)_tors > 18, which excludes quaternionic multiplication (an O-PQM,
//         or non-maximal PQM, genus-2 Jacobian over Q has #tors <= 16, resp. <= 18,
//         by Laga-Schembri-Shnidman-Voight 2024, Thm 1.5 + Rem 1.2.3).
// (ii)+(iii)+(iv) together force End(Jac_Qbar) = Z.  The actual chi polynomials are
// printed so the certificate is checkable without rerunning the search.
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

// returns: is-root-power-strict, chi, splitting-degree
function ChiData(C, p)
  chi := Px!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
  if Degree(chi) ne 4 or not IsIrreducible(chi) then return false, false, chi, 0; end if;
  K<a> := NumberField(chi);
  strict := true;
  for n in [2..12] do
    if Degree(MinimalPolynomial(a^n)) ne 4 then strict := false; break; end if;
  end for;
  return true, strict, chi, Degree(SplittingField(chi));
end function;

curves := [
  <"[2,2,14]",  (x+1)*(x+9)*(2*x-115)*(6*x+55)*(3*x^2-220*x+2652)>,
  <"[2,2,20]",  -(x-1)*(6*x+1)*(2*x+7)*(6217*x+1008)*(21*x^2-161*x+144)>,
  <"[2,4,8]",   x*(3*x-5)*(5*x+333)*(16*x^2+23*x+576)>,
  <"[6,6]",     x*(39*x^2-69*x+125)*(48*x^2+8*x+39)>,
  <"[2,4,4]",   799^2*x*(x+16)*(x+(644/799)^2)*(x^2+x+16)>,
  <"[2,2,4,4]", x*(x+36^2)*(x+57^2)*(x+64^2)*(x+132^2)>,
  <"[2,2,2,8]", x*(x+1)*(x+55^2)*(x+99^2)*(x+125^2)>,
  <"[2,2,2,12]",x*(x-120^2)*(x-143^2)*(x-266^2)*(x-218^2)*(x-241^2)>,
  <"[2,2,2,10]",x*(x+1)*(x-1)*(3*x-7)*(8*x-13)*(24*x+25)>
];

out := "";
for pair in curves do
  lbl := pair[1]; f := pair[2];
  C := HyperellipticCurve(f);
  T := TorsionSubgroup(Jacobian(SimplifiedModel(ReducedMinimalWeierstrassModel(C))));
  inv := Invariants(T); ord := #T;
  D := Integers()!Discriminant(C);
  // (ii) strict root-power prime
  p0 := 0; chi_p := Px!0;
  for p in PrimesInInterval(11, 400) do
    if D mod p eq 0 then continue; end if;
    ok, strict, chi, _ := ChiData(C, p);
    if ok and strict then p0 := p; chi_p := chi; break; end if;
  end for;
  // (iii) disjoint second prime
  _, _, chiP0, dP0 := ChiData(C, p0);
  q0 := 0; chi_q := Px!0;
  for p in PrimesInInterval(11, 400) do
    if p eq p0 or D mod p eq 0 then continue; end if;
    // require root-power strictness at the second prime too: absolute
    // simplicity of the reduction forces End^0 = Q(pi_q), so a real-quadratic
    // center would embed in both splitting fields (rules out even chi like
    // x^4-10x^2+1849, whose reduction splits over the quadratic extension)
    ok, strictq, chi, dQ := ChiData(C, p);
    if not ok or not strictq then continue; end if;
    if Degree(SplittingField(chiP0*chi)) eq dP0*dQ then q0 := p; chi_q := chi; break; end if;
  end for;
  qm := ord gt 18;
  out cat:= Sprintf("CURVE %o\n", lbl);
  out cat:= Sprintf("  torsion = %o (order %o)\n", inv, ord);
  out cat:= Sprintf("  root-power prime p0 = %o : chi = %o  [irreducible, pi^n deg 4 for n<=12 => geometrically simple]\n", p0, chi_p);
  out cat:= Sprintf("  disjoint prime  q0 = %o : chi = %o  [splitting field disjoint from p0 => center(End^0) = Q]\n", q0, chi_q);
  out cat:= Sprintf("  #tors = %o > 18 => not QM [LSSV 2024]  : %o\n", ord, qm);
  out cat:= Sprintf("  ==> End(Jac_Qbar) = Z : %o\n\n", (p0 ne 0) and (q0 ne 0) and qm);
  printf "%o: torsion %o, root-power p0=%o, disjoint q0=%o, tors>18 %o => End=Z %o\n",
    lbl, inv, p0, q0, qm, (p0 ne 0) and (q0 ne 0) and qm;
end for;
Write("data/claude_endz_certificates.txt", out : Overwrite := true);
print "ALL_DONE";
quit;
