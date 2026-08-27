// claude_ov_88modp.m -- Lane 5 (2026-07-25, resumed session).
//
// THE DECISIVE LOCAL TEST for [8,8] on Nicholls' Lambda_334.
//
// If some member C1 of the family over Q had  J1(Q)_tors  containing  Z/8 x Z/8,
// then for every prime p of good reduction the reduction satisfies
//     J1(F_p)  contains  Z/8 x Z/8 .
// The family is the affine 3-fold in (s,t,beta) with  eta^2 = s^2+(1-t^2)(t^2-beta^2)
// and  u = eta/(s*beta)  (verified chart, results/claude_ov_88chart_verify.log,
// TEST1 175/175), and the DOUBLE-STAGE-1 sublocus -- necessary for [8,8], since
// both prescribed 2-torsion classes of J2 must be 2-divisible -- is cut out by the
// three SQUARE conditions (verified against exact Magma torsion, TEST2 235/235):
//     t^2-1 = square,  (t^2-1)(s^2-(t^2-1)^2) = square,  t^2-beta^2 = square.
// Squares reduce to squares, so a global [8,8] member reduces to an F_p-point of
// the same chart satisfying the same conditions with "square" = "QR in F_p".
// Enumerating ALL of (s,t,beta) in F_p^3 is therefore a COMPLETE test of the affine
// chart at p (the only gap is the boundary at infinity of (s,t,beta)-space).
//
// Reported per prime, both over the whole family and over the DS1 sublocus:
//   built     : nonsingular members
//   c44       : J1(F_p) contains (Z/4)^2                (necessary for [8,8])
//   c48       : J1(F_p) contains Z/4 x Z/8
//   c88       : J1(F_p) contains Z/8 x Z/8              <-- the [8,8] gate
// If c88 = 0 on the DS1 sublocus for some p, [8,8] is IMPOSSIBLE on this chart at
// any member with good reduction at p.  The whole-family column is the POSITIVE
// CONTROL: it must be nonzero, otherwise the test is vacuous / miscoded.
//
// Params: PRIMES (comma list), NCH (children), OUTDIR
SetColumns(0);
SetMemoryLimit(4*10^9);
Z := Integers();

if not assigned PRIMES then PRIMES := "11,13,17,19,23,29,31,37,41,43"; end if;
if not assigned NCH then NCH := 10; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned OUTDIR then
  OUTDIR := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/modp";
end if;
System("mkdir -p " cat OUTDIR);
plist := [ StringToInteger(z) : z in Split(PRIMES, ",") ];

function ScanPrime(p)
  F := GF(p);  PF<x> := PolynomialRing(F);
  nbuilt := 0; nds1 := 0;
  c44 := 0; c48 := 0; c88 := 0;
  d44 := 0; d48 := 0; d88 := 0;
  for tv in F do
    if tv eq 0 or tv^2 eq 1 then continue; end if;
    ds1a := IsSquare(tv^2-1);
    for sv in F do
      if sv eq 0 then continue; end if;
      Av := sv^2 - tv^4 + tv^2;
      if Av eq 0 then continue; end if;
      ds1b := ds1a and IsSquare((tv^2-1)*(sv^2-(tv^2-1)^2));
      for bt in F do
        if bt eq 0 then continue; end if;
        e2 := sv^2 + (1-tv^2)*(tv^2 - bt^2);
        if e2 eq 0 then continue; end if;
        isq, eta := IsSquare(e2);
        if not isq then continue; end if;
        uv := eta/(sv*bt);
        if uv^2*sv^2 + 1 - tv^2 eq 0 then continue; end if;
        av := Av/(1 - tv^2);
        bv := bt^2;
        cv := tv^2;
        d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)
                *(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
        if d2 eq 0 then continue; end if;
        l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
        l2 := -x^2 + bv*cv;
        l3 := x^2 - av;
        f1 := l1*l2*l3;
        if Degree(f1) ne 6 then continue; end if;
        ff := d2*f1/LeadingCoefficient(f1);
        if Discriminant(ff) eq 0 then continue; end if;
        nbuilt +:= 1;
        ds1 := ds1b and IsSquare(tv^2 - bt^2);
        if ds1 then nds1 +:= 1; end if;
        G := AbelianGroup(Jacobian(HyperellipticCurve(ff)));
        inv := Invariants(G);
        n4 := #[ i : i in inv | i mod 4 eq 0 ];
        n8 := #[ i : i in inv | i mod 8 eq 0 ];
        if n4 ge 2 then
          c44 +:= 1;  if ds1 then d44 +:= 1; end if;
          if n8 ge 1 then c48 +:= 1; if ds1 then d48 +:= 1; end if; end if;
          if n8 ge 2 then
            c88 +:= 1;  if ds1 then d88 +:= 1; end if;
          end if;
        end if;
      end for;
    end for;
  end for;
  return Sprintf("MODP p=%o built=%o ds1=%o | all: c44=%o c48=%o c88=%o | ds1: d44=%o d48=%o d88=%o",
                 p, nbuilt, nds1, c44, c48, c88, d44, d48, d88);
end function;

printf "MODP_START primes=%o children=%o\n", plist, NCH;
for c in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then
    H := Open(OUTDIR cat "/part" cat IntegerToString(c) cat ".txt", "w");
    i := c+1;
    while i le #plist do
      try
        Puts(H, ScanPrime(plist[i]));
      catch e
        Puts(H, Sprintf("MODP p=%o ERROR %o", plist[i], e`Object));
      end try;
      Flush(H);
      i +:= NCH;
    end while;
    Flush(H); delete H;
    quit;
  end if;
end for;
WaitForAllChildren();
for c in [0..NCH-1] do
  for l in Split(Read(OUTDIR cat "/part" cat IntegerToString(c) cat ".txt"), "\n") do
    if #l gt 0 then printf "%o\n", l; end if;
  end for;
end for;
printf "MODP_DONE\n";
quit;
