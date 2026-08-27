// claude_ov_88lift_profile.m -- Lane 5 (2026-07-25).  [8,8] on Nicholls' Lambda_334.
//
// THE ARGUMENT (this session's structural contribution)
// ----------------------------------------------------
// C1 : y^2 = d2*f1,  f1 = l1*l2*l3 monic sextic.  On the family bc is ALWAYS a square
// (verified 12825/12825), so l2 splits and the factor type is [1,1,2,2]; the rational
// 2-torsion is still (Z/2)^2 = Sigma[2] = {0,[L1],[L2],[L3]} (Galois-stable even
// subsets of orbit type (1,1,2,2) modulo complement).  Hence the 2-primary part of
// J1(Q)_tors has AT MOST TWO cyclic factors, so
//     J1(Q)_tors contains Z/8 x Z/8  ==>  J1(Q)[4] = Sigma  and  Sigma subset 2*J1(Q).
//
// Now: Sigma[2] = 2*Sigma subset 2*J1(Q), so the even-degree x-T homomorphism
//     mu : J1(Q)/2J1(Q) --> L*/Q*(L*)^2,   L = Q[x]/f1,
// is CONSTANT on cosets of Sigma[2].  Sigma/Sigma[2] = (Z/2)^2 and the three nonzero
// cosets are exactly the three classes of the 12 order-4 elements.  Composing mu with
// the component norms gives three characters
//     chi_j : Sigma/Sigma[2] --> Q*/(Q*)^2,   chi_j(D) = Res(l_j^monic, u_D) mod squares
// (well defined: N(lambda*square) = lambda^2*square for the two quadratic components,
// and u(r)u(-r) -> lambda^2 for the split one).  Their common kernel H is a SUBGROUP of
// (Z/2)^2, so #H in {1,2,4} and the number of norm-passing order-4 cosets is 0, 1 or 3.
//     [8,8]  ==>  H = (Z/2)^2, i.e. ALL THREE cosets pass, i.e. chi_j == 1 identically.
// The 2026-07-18 diagnostic saw exactly ONE passing coset (4 of 12 elements).  This
// script measures the character table (chi_j on two coset generators) over the family.
//
// Output per member:
//   CHITAB (s,t,v) | g1: [c1,c2,c3] | g2: [c1,c2,c3] | g3: [c1,c2,c3] | H=<n> | tors
// where c_j are squarefree integers.  H = #passing cosets + 1.
//
// Params: MODE = "box" | "list" | "mn"
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);  Z := Integers();

function SqFree(r)
  assert r ne 0;
  m := Numerator(r)*Denominator(r);
  s := Sign(m);  m := AbsoluteValue(m);
  res := s;
  for pf in Factorization(m) do if IsOdd(pf[2]) then res *:= pf[1]; end if; end for;
  return res;
end function;

function MemberSTV(sv, tv, vv)
  if sv eq 0 or tv eq 0 or tv^2 eq 1 then return false, _, _, _; end if;
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 then return false, _, _, _; end if;
  den := -sv^2*tv*Av*vv^2 + tv;
  if den eq 0 then return false, _, _, _; end if;
  uv := (-sv^2*Av*vv^2 - 2*Av*vv - 1)/den;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false, _, _, _; end if;
  av := Av/(1 - tv^2);
  bv := Av/(uv^2*sv^2 + 1 - tv^2);
  cv := tv^2;
  d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)
          *(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
  if d2 eq 0 then return false, _, _, _; end if;
  l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
  l2 := -x^2 + bv*cv;
  l3 := x^2 - av;
  f1 := l1*l2*l3;
  if Degree(f1) ne 6 then return false, _, _, _; end if;
  f := d2*f1/LeadingCoefficient(f1);
  if Discriminant(f) eq 0 then return false, _, _, _; end if;
  dd := LCM([Denominator(cc) : cc in Coefficients(f)]);
  fint := P![ cc*dd^2 : cc in Coefficients(f) ];
  cont := GCD([Z!cc : cc in Coefficients(fint)]);
  sqf := 1;
  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  fint := P![ (Z!cc) div sqf : cc in Coefficients(fint) ];
  Ls := [ l1/LeadingCoefficient(l1), l2/LeadingCoefficient(l2), l3/LeadingCoefficient(l3) ];
  return true, fint, Ls, [av,bv,cv,Av,uv,d2,sv,tv,vv];
end function;

procedure ChiTable(tag, fint, Ls, aux)
  J := Jacobian(HyperellipticCurve(fint));
  T, mT := TorsionSubgroup(J);
  inv := Invariants(T);
  // collect the character values by coset (keyed by the image of D in T/T[2] ... we
  // simply key by 2*D, which is constant on cosets of Sigma[2] and separates them)
  seen := AssociativeArray();
  bad := 0;
  for g in T do
    if Order(g) ne 4 then continue; end if;
    D := mT(g);  uD := D[1];
    if Degree(uD) ne D[3] then bad +:= 1; continue; end if;
    key := Sprint(mT(2*g)[1]);
    cs := [];
    ok := true;
    for j in [1..3] do
      R := Resultant(Ls[j], uD);
      if R eq 0 then ok := false; break; end if;
      Append(~cs, SqFree(R));
    end for;
    if not ok then bad +:= 1; continue; end if;
    if IsDefined(seen, key) then
      if seen[key] ne cs then printf "INCONSISTENT %o key=%o %o vs %o\n", tag, key, seen[key], cs; end if;
    else
      seen[key] := cs;
    end if;
  end for;
  ks := Sort(Setseq(Keys(seen)));
  npass := #[ k : k in ks | seen[k] eq [1,1,1] ];
  strs := [ Sprintf("%o", seen[k]) : k in ks ];
  printf "CHITAB %o | cosets=%o | %o | passing=%o | tors=%o | bad=%o\n",
         tag, #ks, Join(strs, " | "), npass, inv, bad;
  // auxiliary square classes for pattern fitting
  printf "AUX %o a=%o b=%o c=%o A=%o d2=%o disc1=%o disc3=%o r12=%o r13=%o r23=%o\n",
         tag, SqFree(aux[1]), SqFree(aux[2]), SqFree(aux[3]), SqFree(aux[4]), SqFree(aux[6]),
         SqFree(Discriminant(Ls[1])), SqFree(Discriminant(Ls[3])),
         SqFree(Resultant(Ls[1],Ls[2])), SqFree(Resultant(Ls[1],Ls[3])), SqFree(Resultant(Ls[2],Ls[3]));
end procedure;

if not assigned MODE then MODE := "box"; end if;

if MODE eq "box" then
  if not assigned HB then HB := 3; elif Type(HB) eq MonStgElt then HB := StringToInteger(HB); end if;
  if not assigned SZLIM then SZLIM := 2500; elif Type(SZLIM) eq MonStgElt then SZLIM := StringToInteger(SZLIM); end if;
  if not assigned PART then PART := 0; elif Type(PART) eq MonStgElt then PART := StringToInteger(PART); end if;
  if not assigned NPARTS then NPARTS := 1; elif Type(NPARTS) eq MonStgElt then NPARTS := StringToInteger(NPARTS); end if;
  vals := [Q| ];
  for a in [1..HB] do for b in [1..HB] do
    if GCD(a,b) eq 1 then Include(~vals, Q!a/b); Include(~vals, -Q!a/b); end if;
  end for; end for;
  pos := Sort([ w : w in vals | w gt 0 ]);
  vv2 := Sort(Setseq(Seqset(vals)));
  idx := 0; cnt := 0;
  for sv in pos do for tv in pos do for vw in vv2 do
    idx +:= 1;
    if (idx mod NPARTS) ne PART then continue; end if;
    ok, fint, Ls, aux := MemberSTV(sv, tv, vw);
    if not ok then continue; end if;
    if #Sprint(Coefficients(fint)) gt SZLIM then continue; end if;
    cnt +:= 1;
    tag := Sprintf("(%o,%o,%o)", sv, tv, vw);
    try
      ChiTable(tag, fint, Ls, aux);
    catch e
      printf "ERROR %o %o\n", tag, e`Object;
    end try;
    if cnt mod 25 eq 0 then printf "PROGRESS part=%o count=%o\n", PART, cnt; end if;
  end for; end for; end for;
  printf "BOX_DONE part=%o count=%o\n", PART, cnt;
  quit;
end if;

if MODE eq "mn" then
  // LIST of "m,n,v" in the stage-1 (m,n) coordinates
  if not assigned LIST then printf "need LIST\n"; quit; end if;
  for chunk in Split(LIST, ";") do
    pr := Split(chunk, ",");
    mv := Q!eval(pr[1]); nv := Q!eval(pr[2]); vw := Q!eval(pr[3]);
    al := (mv^2-1)/(2*mv); tv := (mv^2+1)/(2*mv); sv := (nv^2+al^4)/(2*nv);
    ok, fint, Ls, aux := MemberSTV(sv, tv, vw);
    if not ok then printf "INVALID %o\n", chunk; continue; end if;
    tag := Sprintf("mnv=(%o,%o,%o)|stv=(%o,%o,%o)", mv,nv,vw, sv,tv,vw);
    try ChiTable(tag, fint, Ls, aux);
    catch e printf "ERROR %o %o\n", tag, e`Object; end try;
  end for;
  printf "MN_DONE\n";
  quit;
end if;

printf "unknown MODE\n";
quit;
