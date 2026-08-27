// ===========================================================================
// Lane 2 (claude_ov_c7z) -- HOW FAR BEYOND c=2 CAN WE DECIDE?
//
// Slice z_1 = c of the contact-7 order-112 surface S = {e1=6, e5=2, e3=2e4-2}.
// The other four z's are the roots of the pencil
//     P_m(w) = w^4 - e1 w^3 + m w^2 - (A m + B) w + e4,
//     e1 = 6-c, e4 = 2/c, T = 3 - 1/(1-c), A = (T-2)/(T-1),
//     B = (T*(c-5+2/c) - 3c + 14)/(T-1).
// Splitting into two rational quadratics (w^2 - a w + b)(w^2 - (e1-a) w + e4/b)
// forces (a,b) onto
//     Cab_c :  (e1-A-a) b^2 + (A a^2 - A e1 a - B) b + e4 (a - A) = 0,
// a conic in b over the a-line, i.e. the genus-1 curve
//     Cab_c :  y^2 = Delta_c(a) := (A a^2 - A e1 a - B)^2 - 4 (e1-A-a) e4 (a-A).
// (For c=2 Delta has a repeated root and Cab_2 is RATIONAL -- that is exactly the
//  degeneration that made the c=2 slice tractable.)
// Order 112 additionally needs d1 = a^2-4b and d2 = (e1-a)^2-4 e4/b both squares.
//
// STRATEGY: if Cab_c has NO rational point, or genus 1 with RANK 0, the slice is
// DECIDED outright (finitely many (a,b), each tested exactly).  Report which.
// ===========================================================================
SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Pa<aa> := PolynomialRing(Q);

function SliceData(c)
  // returns e1, e4, A, B, Delta (quartic in aa); flag false if degenerate
  if c eq 0 or c eq 1 then return 0,0,0,0,Pa!0,false; end if;
  T := 3 - 1/(1-c);
  if T eq 1 then return 0,0,0,0,Pa!0,false; end if;      // c = 1/2
  e1 := 6-c; e4 := 2/c;
  A := (T-2)/(T-1);
  B := (T*(c-5+2/c) - 3*c + 14)/(T-1);
  D := (A*aa^2 - A*e1*aa - B)^2 - 4*(e1-A-aa)*e4*(aa-A);
  return e1, e4, A, B, D, true;
end function;

// ---- control: c=2 must give a NON-squarefree Delta (Cab_2 rational) ----
e1,e4,A,B,D2,ok := SliceData(2);
printf "CONTROL c=2 : e1=%o e4=%o A=%o B=%o\n", e1, e4, A, B;
printf "CONTROL c=2 : Delta = %o\n", D2;
printf "CONTROL c=2 : Delta factorisation = %o  squarefree=%o\n", Factorisation(D2), IsSquarefree(D2);

clist := [];
for d in [1..4] do
  for n in [-12..12] do
    if n eq 0 then continue; end if;
    if GCD(Abs(n),d) ne 1 then continue; end if;
    cv := Q!n/d;
    if cv in {Q| 0, 1, 1/2 } then continue; end if;
    if not (cv in clist) then Append(~clist, cv); end if;
  end for;
end for;
// the z-values actually occurring in the eleven recorded three-root ([2,2,14]) hits
extra := [ Q| -1/2,4,5/2,-1/9,2,-7/3,-1/4,-8/7,22/7,19/4,9/5,17/21,25/21,18/13,49/39,
            -61/450,625/114,-41/124,-16/17,289/124,297/133,361/525,50/33,189/155,121/155,49/36,50/63 ];
for e in extra do if not (e in clist) then Append(~clist,e); end if; end for;
// -b ONLYN:=2 ONLYD:=3  restricts the run to the single slice c = 2/3
if assigned ONLYN then
  nn := Type(ONLYN) eq MonStgElt select StringToInteger(ONLYN) else ONLYN;
  dd := assigned ONLYD select (Type(ONLYD) eq MonStgElt select StringToInteger(ONLYD) else ONLYD) else 1;
  clist := [ Q!nn/dd ];
end if;
printf "\ndeciding %o slices\n\n", #clist;

ndec := 0; nund := 0; nhit := 0; undec := [];
for c in clist do
  e1,e4,A,B,D,ok := SliceData(c);
  if not ok then printf "c=%-10o : DEGENERATE slice (z-value not allowed)\n", c; continue; end if;
  if D eq 0 then printf "c=%-10o : Delta identically 0\n", c; continue; end if;
  sf := IsSquarefree(D);
  dg := Degree(D);
  if not sf or dg lt 3 then
    printf "c=%-10o : Delta deg %o squarefree=%o  -> Cab_c NOT genus 1 (special slice)\n", c, dg, sf;
    Append(~undec, <c,"special">); nund +:= 1;
    continue;
  end if;
  // integral model of y^2 = Delta: multiply by a SQUARE, so the square class (and
  // hence the set of rational points, up to y -> y/dl) is unchanged
  dl := LCM([ Denominator(x) : x in Coefficients(D) ]);
  Dint := PolynomialRing(Integers()) ! (dl^2 * D);
  C := HyperellipticCurve(Dint);
  g := Genus(C);
  if g ne 1 then
    printf "c=%-10o : genus(Cab_c) = %o (unexpected)\n", c, g;
    Append(~undec, <c,"genus" cat IntegerToString(g)>); nund +:= 1;
    continue;
  end if;
  // (1) is there a rational point at all?
  pts := Points(C : Bound := 3000);
  if #pts eq 0 then
    // no small point: try the local test, which if it fails DECIDES the slice
    lsol := true;
    try lsol := IsLocallySolvable(C); catch e ; end try;
    if not lsol then
      printf "c=%-10o : Cab_c genus 1, NO Qp-point  ==> SLICE DECIDED (empty)\n", c;
      ndec +:= 1;
    else
      printf "c=%-10o : Cab_c genus 1, ELS but no point of height<=3000 -> UNDECIDED\n", c;
      Append(~undec, <c,"nopoint">); nund +:= 1;
    end if;
    continue;
  end if;
  E, mp := EllipticCurve(C, pts[1]);
  E := MinimalModel(E);
  tor := TorsionSubgroup(E);
  lo, hi := RankBounds(E);
  printf "c=%-10o : Cab_c genus 1 -> E %o cond %o  torsion %o  rank %o..%o",
     c, aInvariants(E), Conductor(E), Invariants(tor), lo, hi;
  if hi eq 0 then
    // rank 0: E(Q) = torsion, and C =~ E over Q, so #C(Q) = #E(Q).  Enumerate C(Q)
    // by naive search and CERTIFY completeness by matching the count.
    nE := #TorsionSubgroup(E);
    cpts := Points(C : Bound := 200000);
    printf "\n   rank 0: #E(Q) = %o, naive search on C found %o points -> complete = %o\n",
       nE, #cpts, #cpts eq nE;
    if #cpts ne nE then
      printf "   !! could not exhibit all %o points of C(Q) -- slice NOT decided\n", nE;
      Append(~undec, <c,"rank0-incomplete">); nund +:= 1;
      continue;
    end if;
    nb := 0;
    for pre in cpts do
      if pre[3] eq 0 then continue; end if;
      av := pre[1]/pre[3]; yv := (pre[2]/pre[3]^2)/dl;   // undo the dl^2 scaling
      den := e1 - A - av;
      // den = 0: the equation for b degenerates to LINEAR, one finite root
      // (the other root of the conic has gone to b = infinity, not a configuration)
      bl := [];
      if den eq 0 then
        lin := A*av^2 - A*e1*av - B;
        if lin ne 0 then bl := [ -e4*(av-A)/lin ]; end if;
      else
        bl := [ (-(A*av^2 - A*e1*av - B) + sgn*yv)/(2*den) : sgn in [1,-1] ];
      end if;
      for bv in bl do
        if bv eq 0 then continue; end if;
        // verify it really is on Cab_c
        if den*bv^2 + (A*av^2 - A*e1*av - B)*bv + e4*(av-A) ne 0 then continue; end if;
        nb +:= 1;
        d1 := av^2 - 4*bv;
        d2 := (e1-av)^2 - 4*e4/bv;
        if IsSquare(d1) and IsSquare(d2) then
          printf "\n   *** ORDER-112 CANDIDATE c=%o a=%o b=%o d1=%o d2=%o\n", c, av, bv, d1, d2;
          nhit +:= 1;
        end if;
      end for;
    end for;
    printf "  ==> RANK 0, %o rational (a,b) tested, SLICE DECIDED (no order 112)\n", nb;
    ndec +:= 1;
  else
    printf "  -> rank >= %o, UNDECIDED\n", lo;
    Append(~undec, <c,"rank" cat IntegerToString(lo)>); nund +:= 1;
  end if;
end for;

printf "\n==== SUMMARY ====\n";
printf "slices DECIDED (rank 0 or empty): %o\n", ndec;
printf "slices UNDECIDED               : %o\n", nund;
printf "order-112 candidates found     : %o\n", nhit;
printf "undecided list: %o\n", undec;
printf "CABRANK_DONE\n";
quit;
