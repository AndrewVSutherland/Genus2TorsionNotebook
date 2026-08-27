// ===========================================================================
// Lane 2 (claude_ov_c7z) -- LOCAL point counts on the undecided slice fibres.
//
// For slice z_1 = c the order-112 fibre is
//     X_c :  u^2 = d1 = a^2 - 4b,   v^2 = d2 = (e1-a)^2 - 4 e4/b   over  Cab_c,
// of genus 13 for every c != 2 tested (genus 4 at c=2), and the three-root fibre is
//     D_c :  u^2 = d1                                              over  Cab_c,
// of genus 4 for every c != 2 (genus 1 at c=2).
//
// Weil: a genus-g curve over F_p has a degree-1 place once p+1 > 2 g sqrt(p), i.e.
//   g = 4  -> guaranteed for p >= 59 ;   g = 13 -> guaranteed for p >= 683.
// So an EMPTY count is only possible at small p.  We count degree-1 places of the
// smooth models over F_p for every genus-preserving prime p < 100.  A prime where
// the genus drops is discarded as (a proxy for) bad reduction.
//
// A slice with #X_c(F_p) = 0 at a genus-preserving p is a candidate LOCAL
// OBSTRUCTION -> that slice would be decided outright (modulo checking that the
// model really has good reduction at p, which this script does NOT do).
// ===========================================================================
SetColumns(0);
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();

clist := [Q| 2, -1, 3, 4, 5, -2, 3/2, -1/3, 5/2, 7, 2/3, -5, 10, 6, -1/2, -3, 8, 1/3,
             4/3, -4, 9, -2/3, 11, 12, -6, 5/3, 7/3, 1/4, 3/4, 5/4 ];
plist := [ p : p in [3..97] | IsPrime(p) ];

function SliceCsts(c)
  if c eq 0 or c eq 1 then return false,0,0,0,0; end if;
  T := 3 - 1/(1-c);
  if T eq 1 then return false,0,0,0,0; end if;
  e1 := 6-c; e4 := 2/c;
  A := (T-2)/(T-1);
  B := (T*(c-5+2/c) - 3*c + 14)/(T-1);
  return true, e1, e4, A, B;
end function;

printf "primes tested: %o\n", plist;
printf "%-8o %-6o %-24o %-24o\n", "c", "p", "D_c (3-root, g=4)", "X_c (order-112, g=13)";
hits := [];
for c in clist do
  ok, e1, e4, A, B := SliceCsts(c);
  if not ok then printf "c=%-8o DEGENERATE\n", c; continue; end if;
  // reference genera over Q
  R2<aa,bb> := PolynomialRing(Q,2);
  F := (e1 - A - aa)*bb^2 + (A*aa^2 - A*e1*aa - B)*bb + e4*(aa - A);
  gD0 := -1; gX0 := -1;
  try
    CQ := Curve(AffineSpace(R2), F);
    KQ := FunctionField(CQ); FFQ, mQ := AlgorithmicFunctionField(KQ);
    aQ := mQ(KQ!(R2.1)); bQ := mQ(KQ!(R2.2));
    PQ<uq> := PolynomialRing(FFQ);
    L1Q := FunctionField(uq^2 - (aQ^2 - 4*bQ)); gD0 := Genus(L1Q);
    PQ2<vq> := PolynomialRing(L1Q);
    L2Q := FunctionField(vq^2 - (L1Q!((e1-aQ)^2 - 4*e4/bQ))); gX0 := Genus(L2Q);
  catch e ; end try;
  nz := 0;
  for p in plist do
    dens := [ Denominator(x) : x in [e1,e4,A,B] ];
    if &or[ d mod p eq 0 : d in dens ] then continue; end if;
    k := GF(p);
    Rk<ak,bk> := PolynomialRing(k,2);
    Fk := (k!e1 - k!A - ak)*bk^2 + (k!A*ak^2 - k!A*k!e1*ak - k!B)*bk + k!e4*(ak - k!A);
    try
      Ck := Curve(AffineSpace(Rk), Fk);
      if not IsIrreducible(Fk) then continue; end if;
      Kk := FunctionField(Ck); FFk, mk := AlgorithmicFunctionField(Kk);
      if not IsExact(ConstantField(FFk)) then continue; end if;
      if #ConstantField(FFk) ne p then continue; end if;   // not geometrically irreducible
      aK := mk(Kk!(Rk.1)); bK := mk(Kk!(Rk.2));
      Pk<uk> := PolynomialRing(FFk);
      d1K := aK^2 - 4*bK;
      if IsSquare(d1K) then continue; end if;              // cover degenerates mod p
      L1 := FunctionField(uk^2 - d1K);
      if Genus(L1) ne gD0 then continue; end if;           // genus drop => bad prime
      Pk2<vk> := PolynomialRing(L1);
      d2K := L1!((k!e1-aK)^2 - 4*k!e4/bK);
      if IsSquare(d2K) then continue; end if;
      L2 := FunctionField(vk^2 - d2K);
      if Genus(L2) ne gX0 then continue; end if;
      nD := #Places(L1, 1);
      nX := #Places(L2, 1);
      nz +:= 1;
      if nD eq 0 or nX eq 0 then
        printf "c=%-8o p=%-4o  #D_c(F_p)=%-6o #X_c(F_p)=%-6o   *** EMPTY ***\n", c, p, nD, nX;
        Append(~hits, <c,p,nD,nX>);
      end if;
    catch e ;
    end try;
  end for;
  printf "c=%-8o : genus(D_c)=%o genus(X_c)=%o ; %o genus-preserving primes tested, %o empty\n",
     c, gD0, gX0, nz, #[ h : h in hits | h[1] eq c ];
end for;
printf "\n==== candidate local obstructions ====\n";
if #hits eq 0 then
  printf "NONE: every genus-preserving prime p < 100 on every slice tested has a\n";
  printf "degree-1 place on BOTH D_c and X_c -- no small-prime obstruction exists.\n";
else
  for h in hits do printf "  c=%o p=%o  #D=%o #X=%o\n", h[1], h[2], h[3], h[4]; end for;
end if;
printf "LOCALOBS_DONE\n";
quit;
