// claude_ov_lane7_certscan.m -- Lane 7 (overnight 2026-07-25).
//
// WIDE, FORKED scan for the End=Z certificate primes of all eleven contact-7
// three-root curves.  One forked child per curve; each child sweeps every good
// prime in [11,PMAX], records which are root-power strict, and then searches
// the strict primes for the LEXICOGRAPHICALLY SMALLEST linearly-disjoint pair
// (p0,q0) -- rather than the first pair the serial scan happened to hit -- and
// counts how many disjoint pairs exist among the first NS strict primes
// (a robustness measure: a single disjoint pair is a certificate, many is a
// certificate that cannot be an accident of one prime).
//
// Children write results/claude_ov_lane7_certscan_<i>.txt; the parent
// concatenates them after WaitForAllChildren.
//
// Markers: SCAN / BESTPAIR / LANE7_CERTSCAN_DONE
SetColumns(0);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

if not assigned PMAX then PMAX := 1500; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;
if not assigned NS   then NS := 15;      elif Type(NS)   eq MonStgElt then NS   := StringToInteger(NS);   end if;
if not assigned MemGB then MemGB := 4;   elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;

trips := [
  [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5],
  [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17],
  [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289],
  [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], [-1/2, -13/49, 13/50]
];
tags := ["known", "known", "known-RM", "known", "known", "known",
         "NEW", "NEW", "NEW", "NEW", "NEW"];

outfile := func< i | Sprintf("results/claude_ov_lane7_certscan_%o.txt", i) >;

for i in [1..#trips] do
  pid := Fork();
  if pid ne 0 then continue; end if;   // parent: next child

  // ---------------- child ----------------
  SetMemoryLimit(MemGB*10^9);
  s := trips[i][1]; t := trips[i][2]; u := trips[i][3];
  c4 := (Gfun(s)-Gfun(t))/(s^2-t^2); c0 := Gfun(s) - c4*s^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  f := (h^2 + (x-1)^7) div x^2;
  den := LCM([Denominator(co) : co in Coefficients(f)]);
  F := P![ co*den^2 : co in Coefficients(f) ];
  C := HyperellipticCurve(F);
  D0 := Z!Discriminant(C);
  K2 := QuadraticField(2);

  gp := []; gchi := []; gsf := []; nirr := 0; ngood := 0; nsqrt2 := 0;
  for p in PrimesInInterval(11, PMAX) do
    if D0 mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 then continue; end if;
    ngood +:= 1;
    if not IsIrreducible(chi) then continue; end if;
    nirr +:= 1;
    KK<aa> := NumberField(chi); st := true;
    for n in [2..12] do
      if Degree(MinimalPolynomial(aa^n)) ne 4 then st := false; break; end if;
    end for;
    if not st then continue; end if;
    SF := SplittingField(chi);
    if IsSubfield(K2, SF) then nsqrt2 +:= 1; end if;
    Append(~gp, p); Append(~gchi, chi); Append(~gsf, Degree(SF));
  end for;

  // best (smallest) disjoint pair among the strict primes, + count over first NS
  ns := Min(NS, #gp);
  bestp := 0; bestq := 0; bchi := P!0; bchq := P!0; ndisj := 0;
  for ii in [1..ns] do
    for jj in [ii+1..ns] do
      if Degree(SplittingField(gchi[ii]*gchi[jj])) eq gsf[ii]*gsf[jj] then
        ndisj +:= 1;
        if bestp eq 0 then
          bestp := gp[ii]; bestq := gp[jj]; bchi := gchi[ii]; bchq := gchi[jj];
        end if;
      end if;
    end for;
  end for;

  out := Sprintf("SCAN %o (%o, %o, %o) PMAX=%o goodprimes=%o irreducible=%o rootpower_strict=%o strict_with_sqrt2=%o\n",
                 tags[i], s, t, u, PMAX, ngood, nirr, #gp, nsqrt2);
  out cat:= Sprintf("  strict primes (first 40): %o\n", gp[1..Min(40,#gp)]);
  out cat:= Sprintf("BESTPAIR (%o, %o, %o) p0=%o q0=%o  disjoint_pairs_among_first_%o_strict = %o of %o\n",
                    s, t, u, bestp, bestq, ns, ndisj, (ns*(ns-1)) div 2);
  out cat:= Sprintf("  chi_p0 = %o\n  chi_q0 = %o\n", bchi, bchq);
  Write(outfile(i), out : Overwrite := true);
  quit;
end for;

WaitForAllChildren();
for i in [1..#trips] do
  fh := Open(outfile(i), "r");
  txt := "";
  repeat
    l := Gets(fh);
    if not IsEof(l) then txt cat:= l cat "\n"; end if;
  until IsEof(l);
  printf "%o", txt;
end for;
print "LANE7_CERTSCAN_DONE";
quit;
