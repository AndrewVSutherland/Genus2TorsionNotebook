// claude_2rank_odd_search.m — direct funnel for GENERIC [2,2,14] and [2,22].
// [2,2,14]: f = x(x-1)(x-a)(x-b)(x^2+cx+d), factor type [1,1,1,1,2] => 2-rank 3;
//           hard gate 7 | #J(F_p); exact TorsionSubgroup; End=Z certificate.
// [2,22]:   f = x(x-1)(x^2+ax+b)(x^2+cx+d), factor type [1,1,2,2] => 2-rank 2;
//           hard gate 11 | #J(F_p).
// Params: TARGET ("2214"|"222"), H (box), NP (gate primes), P0,PN (shard on a).
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned TARGET then TARGET := "2214"; end if;
if not assigned H then H := 24; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NP then NP := 7; elif Type(NP) eq MonStgElt then NP := StringToInteger(NP); end if;
if not assigned P0 then P0 := -H; elif Type(P0) eq MonStgElt then P0 := StringToInteger(P0); end if;
if not assigned PN then PN := H; elif Type(PN) eq MonStgElt then PN := StringToInteger(PN); end if;
Q := Rationals(); P<x> := PolynomialRing(Q);

ell := TARGET eq "2214" select 7 else 11;
wantinv := TARGET eq "2214" select [2,2,14] else [2,22];
gateprimes := [ p : p in PrimesInInterval(13, 500) | p ne ell ];

function EndZ(C, ord)
  D := Integers()!Discriminant(C); p0 := 0; chip := P!0; dP := 0; nf := 0;
  for p in PrimesInInterval(11, 300) do
    if D mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<a> := NumberField(chi); st := true;
    for n in [2..12] do if Degree(MinimalPolynomial(a^n)) ne 4 then st := false; break; end if; end for;
    dL := Degree(SplittingField(chi));
    if st and p0 eq 0 then p0 := p; chip := chi; dP := dL; continue; end if;
    if p0 ne 0 and Degree(SplittingField(chip*chi)) eq dP*dL then return (ord gt 18), p0, p, chip, chi; end if;
    nf +:= 1; if nf ge 14 then break; end if;
  end for;
  return false, 0, 0, P!0, P!0;
end function;

tested := 0; gated := 0; hits := 0;
t0 := Cputime();
if TARGET eq "2214" then
  for a in [P0..PN] do
   if a in [0,1] then continue; end if;
   for b in [-H..H] do
    if b in [0,1,a] then continue; end if;
    for c in [-H..H] do for d in [-H..H] do
      if c^2-4*d ge 0 and IsSquare(c^2-4*d) then continue; end if;  // keep quadratic irreducible
      if d eq 0 then continue; end if;
      f := x*(x-1)*(x-a)*(x-b)*(x^2+c*x+d);
      if Degree(f) ne 6 then continue; end if;
      dsc := Discriminant(f); if dsc eq 0 then continue; end if;
      tested +:= 1;
      good := 0; ok := true;
      for p in gateprimes do
        if good ge NP then break; end if;
        if Integers()!dsc mod p eq 0 then continue; end if;
        np := #Jacobian(HyperellipticCurve(PolynomialRing(GF(p))!f));
        good +:= 1; if np mod ell ne 0 then ok := false; break; end if;
      end for;
      if not ok or good lt NP then continue; end if;
      gated +:= 1;
      T := TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
      inv := Invariants(T);
      printf "GATE a=%o b=%o c=%o d=%o tors=%o\n", a,b,c,d, inv;
      if inv eq wantinv then
        ez, p0, q0, cp, cq := EndZ(HyperellipticCurve(f), &*wantinv);
        printf "HIT %o a=%o b=%o c=%o d=%o End=Z:%o\n f=%o\n", TARGET, a,b,c,d, ez, f;
        if ez then printf "GENERIC_HIT_%o f=%o chi_p0(%o)=%o chi_q0(%o)=%o\n", TARGET, f, p0,cp,q0,cq; hits +:= 1; end if;
      end if;
    end for; end for;
   end for;
   printf "PROGRESS a=%o tested=%o gated=%o hits=%o t=%o\n", a, tested, gated, hits, Cputime(t0);
  end for;
else
  for a in [P0..PN] do
   for b in [-H..H] do
    if b eq 0 or (a^2-4*b ge 0 and IsSquare(a^2-4*b)) then continue; end if;
    for c in [-H..H] do for d in [-H..H] do
      if d eq 0 or (c^2-4*d ge 0 and IsSquare(c^2-4*d)) then continue; end if;
      f := x*(x-1)*(x^2+a*x+b)*(x^2+c*x+d);
      if Degree(f) ne 6 then continue; end if;
      dsc := Discriminant(f); if dsc eq 0 then continue; end if;
      tested +:= 1;
      good := 0; ok := true;
      for p in gateprimes do
        if good ge NP then break; end if;
        if Integers()!dsc mod p eq 0 then continue; end if;
        np := #Jacobian(HyperellipticCurve(PolynomialRing(GF(p))!f));
        good +:= 1; if np mod ell ne 0 then ok := false; break; end if;
      end for;
      if not ok or good lt NP then continue; end if;
      gated +:= 1;
      inv := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(f))));
      printf "GATE a=%o b=%o c=%o d=%o tors=%o\n", a,b,c,d, inv;
      if inv eq wantinv then
        ez, p0, q0, cp, cq := EndZ(HyperellipticCurve(f), &*wantinv);
        printf "HIT %o a=%o b=%o c=%o d=%o End=Z:%o\n f=%o\n", TARGET, a,b,c,d, ez, f;
        if ez then printf "GENERIC_HIT_%o f=%o chi_p0(%o)=%o chi_q0(%o)=%o\n", TARGET, f, p0,cp,q0,cq; hits +:= 1; end if;
      end if;
    end for; end for;
   end for;
   printf "PROGRESS a=%o tested=%o gated=%o hits=%o t=%o\n", a, tested, gated, hits, Cputime(t0);
  end for;
end if;
printf "SEARCH_DONE TARGET=%o tested=%o gated=%o hits=%o t=%o\n", TARGET, tested, gated, hits, Cputime(t0);
quit;
