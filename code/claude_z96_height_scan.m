// claude_z96_height_scan.m — direct hunt for Z/96 = 32*3 members of the reconstructed
// Elkies [32] family: enumerate parameters u of bounded height, gate by the necessary
// condition 3 | #J(F_p) at several good primes (finite-prefilters), then exact 3-torsion
// (TorsionSubgroup) on survivors.  Complements the symbolic contact-cover analysis.
// Params: HB (u-height bound, default 120), NP (gate primes, default 8)
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned HB then HB := 120; elif Type(HB) eq MonStgElt then HB := StringToInteger(HB); end if;
if not assigned NP then NP := 8; elif Type(NP) eq MonStgElt then NP := StringToInteger(NP); end if;
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
Q := Rationals(); Pq<t> := PolynomialRing(Q);

cfs := [ Numerator(Coefficient(ftil, i)) : i in [0..5] ];
gateprimes := [ p : p in PrimesInInterval(11, 200) ];

function MemberPoly(uv)
  return Pq![ Evaluate(c, uv) : c in cfs ];
end function;

// height enumeration of u
tested := 0; passed := 0; hits := 0;
for hgt in [1..HB] do
  vals := [];
  for num in [-hgt..hgt] do
    if GCD(num, hgt) ne 1 then continue; end if;
    Append(~vals, Q!num/hgt);
    if num ne 0 and hgt gt 1 then Append(~vals, Q!hgt/num); end if;
  end for;
  for uv in vals do
    fv := MemberPoly(uv);
    if Degree(fv) ne 5 or Discriminant(fv) eq 0 then continue; end if;
    tested +:= 1;
    // gate: 3 | #J(F_p) for NP good primes
    den := LCM([Denominator(c) : c in Coefficients(fv)]);
    good := 0; ok := true;
    for p in gateprimes do
      if good ge NP then break; end if;
      if den mod p eq 0 then continue; end if;
      fp := PolynomialRing(GF(p))![ GF(p) | Numerator(c)*InverseMod(Denominator(c), p) mod p : c in Coefficients(fv) ];
      if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
      np := #Jacobian(HyperellipticCurve(fp));
      good +:= 1;
      if np mod 3 ne 0 then ok := false; break; end if;
    end for;
    if not ok or good lt NP then
      if ok and good lt NP then printf "LOWGOOD u=%o good=%o\n", uv, good; end if;
      continue;
    end if;
    passed +:= 1;
    printf "GATEPASS u=%o (testing exact torsion)\n", uv;
    fexact := Pq![ c*den^2 : c in Coefficients(fv) ];   // integral model (y scaled by den)
    T := TorsionSubgroup(Jacobian(HyperellipticCurve(fexact)));
    inv := Invariants(T);
    printf "EXACT u=%o torsion=%o\n", uv, inv;
    if #T mod 3 eq 0 then
      hits +:= 1;
      printf "HIT_Z96 u=%o torsion=%o f=%o\n", uv, inv, fv;
    end if;
  end for;
  if hgt mod 20 eq 0 then printf "PROGRESS height=%o tested=%o passed=%o hits=%o\n", hgt, tested, passed, hits; end if;
end for;
printf "SEARCH_DONE HB=%o tested=%o gatepassed=%o hits=%o\n", HB, tested, passed, hits;
quit;
