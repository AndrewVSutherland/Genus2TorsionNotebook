// Lane 8 (overnight 2026-07-25), part C: full verification of the CANONICAL
// smallest model recovered for each HLP split target, plus the order-64 bonus
// pair E_{2,8} x E_{2,8}  ([8,8] and [2,4,8], both open targets).
//
// For each model: exact TorsionSubgroup, 2-torsion, factor type, G2 invariants,
// the two explicit degree-2 elliptic subcovers (u = x^2 and u = 1/x^2, available
// because every Prop-4 model is even in x), and the geometric-simplicity /
// RM screen over ALL good primes below PMAX (forked over children so the whole
// prime range is scanned, not just the first success).
SetColumns(0);
if not assigned NCH  then NCH := 12;   elif Type(NCH)  eq MonStgElt then NCH  := StringToInteger(NCH);  end if;
if not assigned PMAX then PMAX := 500; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

QQ := Rationals();
P<x> := PolynomialRing(QQ);
Rz<z> := PolynomialRing(QQ);

models := [
  <"Z/63      HLP Forum Math 12 (2000) eq.(4)  [TRANSCRIBED, verified]",
   897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496>,
  <"Z/7xZ/7   HLP Forum Math 12 (2000) p.15    [TRANSCRIBED, verified]",
   x^6 + 3025*x^4 + 3232987*x^2 + 869675859>,
  <"Z/35      E_7^{-1} glued to E_5^{13/2}     [DERIVED here via Prop 4]",
   9295*x^6 + 2480*x^4 + 5040*x^2 + 640>,
  <"Z/45      E_9^{-5} glued to E_5^{93/10}    [DERIVED here via Prop 4]",
   13981*x^6 + 29240200*x^4 + 49996210000*x^2 + 168300000000>,
  <"Z/5xZ/10  E_10^{-1} glued to E_10^{4/5}    [DERIVED here via Prop 4]",
   45600*x^6 + 289161*x^4 + 35186670*x^2 - 705688215>,
  <"Z/2xZ/24  E_{2,6}^{-1} glued to E_{2,8}^{1/4} [DERIVED here via Prop 4]",
   840*x^6 + 409*x^4 + 46*x^2 + 1>,
  <"Z/3xZ/12  E_12^{1/3} glued to E_6^{-1/4}   [DERIVED here via Prop 4]",
   1296*x^6 - 936*x^4 + 249*x^2 - 20>,
  <"Z/70      Howe's f70 (integral model 144*f70), quintic [TRANSCRIBED, verified]",
   3168*x^5 + 697*x^4 - 23220*x^3 + 37620*x^2 - 23328*x + 5184>
];

dir := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/lane8c";
System("rm -rf " cat dir cat "; mkdir -p " cat dir);

for c in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then
    Fh := Open(dir cat "/part" cat IntegerToString(c) cat ".txt", "w");
    for i in [c+1..#models by NCH] do
      nm := models[i][1]; f := models[i][2];
      Puts(Fh, "==================================================");
      Puts(Fh, nm);
      Puts(Fh, "MODEL  y^2 = " cat Sprint(f));
      C := HyperellipticCurve(f);
      J := Jacobian(C);
      T := TorsionSubgroup(J);
      Puts(Fh, "  genus            : " cat Sprint(Genus(C)));
      Puts(Fh, "  TORSION          : " cat Sprint(Invariants(T)) cat "   order " cat IntegerToString(#T));
      Puts(Fh, "  2-torsion        : " cat Sprint(Invariants(TwoTorsionSubgroup(J))));
      fac := Factorization(f);
      unit := f div &*[t[1]^t[2] : t in fac];
      assert unit*&*[t[1]^t[2] : t in fac] eq f;
      Puts(Fh, "  factor type      : " cat Sprint([Degree(t[1]) : t in fac]) cat "  (leading unit " cat Sprint(unit) cat ")");
      dsc := Integers()!Discriminant(f);
      Puts(Fh, "  disc(f) factored : " cat Sprint(Factorization(dsc)));
      Puts(Fh, "  G2Invariants     : " cat Sprint(G2Invariants(C)));
      Cr, mr := ReducedModel(SimplifiedModel(C));
      Puts(Fh, "  reduced model    : " cat Sprint(Cr));
      // explicit degree-2 elliptic subcovers when f is even in x
      cs := Coefficients(f);
      if Degree(f) eq 6 and &and[cs[k] eq 0 : k in [2..#cs by 2]] then
        cc := [cs[1],cs[3],cs[5],cs[7]];
        MkE := function(A,B,Cq,Dq) return EllipticCurve(z^3 + B*z^2 + A*Cq*z + A^2*Dq); end function;
        E1 := MkE(cc[4],cc[3],cc[2],cc[1]);
        E2 := MkE(cc[1],cc[2],cc[3],cc[4]);
        for e in [<"E1 (u = x^2)", E1>, <"E2 (u = 1/x^2)", E2>] do
          Puts(Fh, "  SUBCOVER " cat e[1] cat " : " cat Sprint(aInvariants(e[2]))
                   cat "  cond " cat IntegerToString(Conductor(e[2]))
                   cat "  torsion " cat Sprint(Invariants(TorsionSubgroup(e[2])))
                   cat "  rank " cat IntegerToString(Rank(e[2])));
        end for;
        okp := true; nchk := 0;
        for p in PrimesInInterval(3, PMAX) do
          if dsc mod p eq 0 then continue; end if;
          Cp := ChangeRing(C, GF(p));
          if not IsNonsingular(Cp) then continue; end if;
          if Evaluate(LPolynomial(Cp),1) ne #ChangeRing(E1,GF(p))*#ChangeRing(E2,GF(p)) then okp := false; end if;
          nchk +:= 1;
        end for;
        Puts(Fh, "  SPLIT CHECK  #J(F_p) = #E1(F_p)*#E2(F_p) at " cat IntegerToString(nchk) cat " good primes: " cat Sprint(okp));
      end if;
      // geometric-simplicity / RM screen over the WHOLE prime range
      nirr := 0; ntest := 0; cores := {}; ftypes := {};
      for p in PrimesInInterval(3, PMAX) do
        if dsc mod p eq 0 then continue; end if;
        Cp := ChangeRing(C, GF(p));
        if not IsNonsingular(Cp) then continue; end if;
        chi := Rz!Reverse(Coefficients(LPolynomial(Cp)));
        ntest +:= 1;
        if IsIrreducible(chi) then nirr +:= 1; end if;
        Include(~ftypes, Sort([Degree(t[1]) : t in Factorization(chi)]));
        dd := Integers()!(Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p));
        Include(~cores, dd eq 0 select 0 else Squarefree(dd));
      end for;
      Puts(Fh, "  good primes < " cat IntegerToString(PMAX) cat " : " cat IntegerToString(ntest)
               cat "   IRREDUCIBLE chi count: " cat IntegerToString(nirr));
      Puts(Fh, "  chi factor-degree types seen : " cat Sprint(Sort(Setseq(ftypes))));
      Puts(Fh, "  real-subfield squarefree cores (RM screen) : " cat Sprint(Sort(Setseq(cores))));
      Puts(Fh, "  VERDICT: " cat (nirr eq 0 select "NOT geometrically simple (chi reducible at every good prime tested) -- SPLIT, as expected"
                                              else "chi irreducible at some prime -- needs the strict two-prime certificate"));
    end for;
    Flush(Fh); delete Fh;
    quit;
  end if;
end for;
WaitForAllChildren();

for c in [0..NCH-1] do
  fn := dir cat "/part" cat IntegerToString(c) cat ".txt";
  txt := Read(fn);
  for l in Split(txt, "\n") do if #l gt 0 then print l; end if; end for;
end for;
print "SEARCH_DONE verify";
quit;
