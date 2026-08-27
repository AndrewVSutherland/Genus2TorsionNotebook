// claude_ov_88lift_worker.m -- Lane 5, sharded worker for the [8,8] lift-layer profile
// on Nicholls' Lambda_334 family.  Reads a member list (one "s t v" or "m n v" per
// line, first token = "stv" or "mnv") and processes the lines with
// (index mod NPARTS) eq PART.
//
// For each member it computes:
//   * the irreducible factorization type of the sextic f1 (always [1,1,2,2] generically)
//   * TorsionSubgroup(J1) exactly
//   * for each of the 12 order-4 elements D of Sigma:
//       - which rational 2-torsion class T = 2D is (index into [L1,L2,L3])
//       - the norm presieve pattern: is Res(L_j^monic, u_D) a square, j=1,2,3
//       - the squarefree class of each Res  (the OBSTRUCTION class when nonsquare)
//       - when the presieve passes: the per-component rational representatives
//         lambda_i (Hilbert 90) and whether a common lambda exists (= full x-T test)
//
// Params: MEMFILE, PART, NPARTS, TORSTIMEOUT (seconds, default 600)
SetColumns(0);
SetMemoryLimit(2*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);
Z := Integers();

if not assigned MEMFILE then printf "need MEMFILE\n"; quit; end if;
if not assigned PART then PART := 0; elif Type(PART) eq MonStgElt then PART := StringToInteger(PART); end if;
if not assigned NPARTS then NPARTS := 1; elif Type(NPARTS) eq MonStgElt then NPARTS := StringToInteger(NPARTS); end if;

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

function SqFree(r)
  if r eq 0 then return 0; end if;
  m := Numerator(r)*Denominator(r);
  s := Sign(m);  m := AbsoluteValue(m);
  res := s;
  for pf in Factorization(m) do
    if IsOdd(pf[2]) then res *:= pf[1]; end if;
  end for;
  return res;
end function;

// full x-T test given the irreducible components
function xTFull(uD, facs)
  lams := []; ds := [];
  for g in facs do
    if Degree(g) eq 1 then
      r := -Coefficient(g,0);
      gam := Evaluate(uD, r);
      if gam eq 0 then return false, 0, [], []; end if;
      Append(~lams, SqFree(gam));  Append(~ds, 1);
    else
      K := NumberField(g);  th := K.1;
      gam := Evaluate(uD, th);
      if gam eq 0 then return false, 0, [], []; end if;
      isq, ng := IsSquare(Norm(gam));
      if not isq then return false, 0, [], []; end if;
      del := 1 + gam/ng;
      if del eq 0 then ng := -ng; del := 1 + gam/ng; end if;
      lam := ng/Norm(del);
      Append(~lams, SqFree(lam));
      Append(~ds, SqFree(Discriminant(g)));
    end if;
  end for;
  k := #lams;
  for mask in [0 .. 2^k - 1] do
    cs := [ SqFree(lams[i]*(((mask div 2^(i-1)) mod 2) eq 0 select 1 else ds[i])) : i in [1..k] ];
    if #Seqset(cs) eq 1 then return true, cs[1], lams, ds; end if;
  end for;
  return false, 0, lams, ds;
end function;

lines := Split(Read(MEMFILE), "\n");
idx := -1;
for ln in lines do
  if ln eq "" or ln[1] eq "#" then continue; end if;
  idx +:= 1;
  if idx mod NPARTS ne PART then continue; end if;
  tk := Split(ln, " ");
  kind := tk[1];
  p1 := Q!eval(tk[2]); p2 := Q!eval(tk[3]); p3 := Q!eval(tk[4]);
  if kind eq "mnv" then
    al := (p1^2-1)/(2*p1); tv := (p1^2+1)/(2*p1);
    if p2 eq 0 then printf "MEMBER %o INVALID\n", ln; continue; end if;
    sv := (p2^2 + al^4)/(2*p2); vv := p3;
  else
    sv := p1; tv := p2; vv := p3;
  end if;
  ok, fint, Ls, aux := MemberSTV(sv, tv, vv);
  if not ok then printf "MEMBER %o INVALID\n", ln; continue; end if;
  if #Sprint(Coefficients(fint)) gt 6000 then printf "MEMBER %o TOOBIG\n", ln; continue; end if;
  try
    C := HyperellipticCurve(fint);  J := Jacobian(C);
    facs := [ g[1]/LeadingCoefficient(g[1]) : g in Factorization(fint) ];
    shape := Sort([ Degree(g) : g in facs ]);
    T, mT := TorsionSubgroup(J);
    inv := Invariants(T);
    // the three rational 2-torsion classes as Jacobian points, in the L1,L2,L3 order
    Tpts := [];
    for j in [1..3] do
      L := Ls[j];
      cofac := fint div (LeadingCoefficient(fint) * &*[Ls[i] : i in [1..3]]) ;  // constant
      Append(~Tpts, L);
    end for;
    n4 := 0; nnorm := 0; npass := 0;
    recs := [];
    for g in T do
      if Order(g) ne 4 then continue; end if;
      n4 +:= 1;
      D := mT(g);  D2 := mT(2*g);
      uD := D[1];
      if Degree(uD) ne D[3] or Degree(uD) ne 2 then
        Append(~recs, <-1, "XXX", [Z|0,0,0], false>);  continue;
      end if;
      // which T_j is 2D?
      wj := 0;
      u2 := D2[1]/LeadingCoefficient(D2[1]);
      for j in [1..3] do if u2 eq Ls[j] then wj := j; end if; end for;
      res := [Z| ];  nsq := [];
      degen := false;
      for j in [1..3] do
        R := Resultant(Ls[j], uD);
        if R eq 0 then degen := true; break; end if;
        Append(~res, SqFree(R));  Append(~nsq, IsSquare(R));
      end for;
      if degen then Append(~recs, <wj, "DEG", [Z|0,0,0], false>); continue; end if;
      pat := &cat[ (b select "1" else "0") : b in nsq ];
      cp := false;
      if pat eq "111" then
        nnorm +:= 1;
        cp, wit, lams, ds := xTFull(uD, facs);
        if cp then npass +:= 1; end if;
      end if;
      Append(~recs, <wj, pat, res, cp>);
    end for;
    printf "MEMBER %o shape=%o tors=%o order4=%o norm111=%o xTpass=%o\n",
           ln, shape, inv, n4, nnorm, npass;
    for r in recs do
      printf "  D4 above=T%o pat=%o resclass=%o xT=%o\n", r[1], r[2], r[3], r[4];
    end for;
  catch e
    printf "MEMBER %o ERROR %o\n", ln, e`Object;
  end try;
end for;
printf "PART_DONE part=%o\n", PART;
quit;
