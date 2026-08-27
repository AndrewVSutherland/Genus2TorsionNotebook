// claude_ov_88gate_validate.m -- Lane 5 (2026-07-25, resumed session).
//
// VALIDATE the [8,8]-containment gate used by the whole [8,8] lift layer.
//
// The gate.  For a genus-2 curve C: y^2 = F(x) over Q with L = Q[x]/F, the even-degree
// x-T homomorphism  mu : J(Q)/2J(Q) --> L^*/(L^*)^2 Q^*  is INJECTIVE.  For an element
// D of J(Q) of order 4 we therefore have
//     D in 2 J(Q)   <=>   mu(D) = 1   <=>   exists lambda in Q^* with
//                          u_D(theta_i) in lambda*(K_i^*)^2 for every component K_i.
// And D = 2E with D of order 4 forces E of order exactly 8.  So
//     "gate passes for D"  <=>  "D has a rational half, i.e. Z/8 sits above D".
// If the rational 2-torsion has rank 2 (factor type [1,1,2,2] on our family), the
// 2-primary part of J(Q)_tors has at most two cyclic factors, hence
//     J(Q)_tors contains Z/8 x Z/8  <=>  ALL 12 order-4 elements of J(Q)[4]=(Z/4)^2
//                                        pass the gate.
// A [4,8] curve must give EXACTLY 4 (one coset of J[2]); a [4,4] curve exactly 0.
// THAT is the discrimination the overnight run never checked: all 32 recorded
// gate-passes were on [4,4] curves, so the gate had never been seen to return
// anything but 0 on a real curve.
//
// GROUND TRUTH.  Because any rational half of a torsion point is itself torsion,
// "D in 2 J(Q)" for D torsion is decided EXACTLY by TorsionSubgroup: D in 2*T.
// This script compares gate vs ground truth on LMFDB curves whose torsion is
// [4,4] (expect 0/12), [4,8], [4,16], [2,4,8], [4,12], [2,2,16], [3,24], [2,48].
//
// Params: FILE (default data/claude_ov_88gate_curves.txt)
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);  Z := Integers();

if not assigned FILE then FILE := "data/claude_ov_88gate_curves.txt"; end if;

function SqFree(r)
  if r eq 0 then return 0; end if;
  m := Numerator(r)*Denominator(r);
  s := Sign(m);  m := AbsoluteValue(m);
  res := s;
  for pf in Factorization(m) do if IsOdd(pf[2]) then res *:= pf[1]; end if; end for;
  return res;
end function;

// The GATE, exactly as used by claude_ov_88lift_worker.m / claude_ov_88conic_base.m,
// but keyed off the irreducible factorization (no assumption of three quadratics).
// returns  status, lams, ds, compat
function Gate(uD, facs)
  lams := []; ds := [];
  for g in facs do
    if Degree(g) eq 1 then
      r := -Coefficient(g,0);
      gam := Evaluate(uD, r);
      if gam eq 0 then return "degenerate", [], [], false; end if;
      Append(~lams, SqFree(gam));  Append(~ds, 1);
    else
      K := NumberField(g);  th := K.1;
      gam := Evaluate(uD, th);
      if gam eq 0 then return "degenerate", [], [], false; end if;
      isq, ng := IsSquare(Norm(gam));
      if not isq then return "normfail", [], [], false; end if;   // presieve failure
      del := 1 + gam/ng;
      if del eq 0 then ng := -ng; del := 1 + gam/ng; end if;
      lam := ng/Norm(del);
      if gam ne lam*del^2 then return "h90fail", [], [], false; end if;
      Append(~lams, SqFree(lam));  Append(~ds, SqFree(Discriminant(g)));
    end if;
  end for;
  k := #lams;
  for mask in [0 .. 2^k - 1] do
    cs := [ SqFree(lams[i]*(((mask div 2^(i-1)) mod 2) eq 0 select 1 else ds[i])) : i in [1..k] ];
    if #Seqset(cs) eq 1 then return "ok", lams, ds, true; end if;
  end for;
  return "ok", lams, ds, false;
end function;

lines := Split(Read(FILE), "\n");
nc := 0; nagree := 0; ndis := 0; ntot4 := 0;
for ln in lines do
  if ln eq "" or ln[1] eq "#" then continue; end if;
  pr := Split(ln, "|");
  lab := pr[1];  while lab[#lab] eq " " do lab := Substring(lab,1,#lab-1); end while;
  fc := [ Q!StringToInteger(z) : z in Split(pr[2], ",") ];
  hc := [ Q!StringToInteger(z) : z in Split(pr[3], ",") ];
  f := P!fc;  h := P!hc;
  F := 4*f + h^2;
  // strip square content
  dd := LCM([Denominator(cc) : cc in Coefficients(F)]);
  F := P![ Z!(cc*dd^2) : cc in Coefficients(F) ];
  cont := GCD([Z!cc : cc in Coefficients(F)]);
  sqf := 1;  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  F := P![ (Z!cc) div sqf : cc in Coefficients(F) ];
  nc +:= 1;
  try
    C := HyperellipticCurve(F);  J := Jacobian(C);
    T, mT := TorsionSubgroup(J);
    inv := Invariants(T);
    facs := [ g[1]/LeadingCoefficient(g[1]) : g in Factorization(F) ];
    shape := Sort([ Degree(g) : g in facs ]);
    // the image 2*T, as a set of Jacobian points
    dbl := { mT(2*g) : g in T };
    n4 := 0; gpass := 0; truth := 0; ndegu := 0; nnormfail := 0; bad := 0;
    nhand := 0; hgate := 0; htruth := 0;
    for g in T do
      if Order(g) ne 4 then continue; end if;
      n4 +:= 1;
      D := mT(g);  uD := D[1];
      tr := D in dbl;
      if tr then truth +:= 1; end if;
      // The gate is only defined on classes whose Mumford u has full degree; the
      // family code (claude_ov_88lift_worker / conic_base) rejects the others, and
      // they never occur on Lambda_334.  Count them separately here.
      if Degree(uD) ne D[3] then ndegu +:= 1; end if;
      if Degree(uD) eq D[3] then nhand +:= 1; if tr then htruth +:= 1; end if; end if;
      st, lams, ds, cp := Gate(uD, facs);
      if Degree(uD) eq D[3] and cp then hgate +:= 1; end if;
      if st eq "normfail" then nnormfail +:= 1; end if;
      if st eq "degenerate" then bad +:= 1; continue; end if;
      if cp then gpass +:= 1; end if;
      if cp ne tr then
        printf "  DISAGREE %o deg(F)=%o u=%o degu=%o d=%o gate=%o truth=%o lam=%o ds=%o\n",
               lab, Degree(F), uD, Degree(uD), D[3], cp, tr, lams, ds;
        if Degree(uD) eq D[3] then ndis +:= 1; end if;
      end if;
    end for;
    ntot4 +:= nhand;
    ag := (hgate eq htruth);
    if ag then nagree +:= 1; end if;
    printf "GATE %o degF=%o shape=%o tors=%o order4=%o handled=%o gate_pass=%o truth_pass=%o normfail=%o degu_mismatch=%o AGREE=%o\n",
           lab, Degree(F), shape, inv, n4, nhand, hgate, htruth, nnormfail, ndegu, ag;
  catch e
    printf "GATE %o ERROR %o\n", lab, e`Object;
  end try;
end for;
printf "GATE_SUMMARY curves=%o order4_classes=%o per_curve_agree=%o disagreements=%o\n",
       nc, ntot4, nagree, ndis;
printf "GATE_VALIDATE_DONE\n";
quit;
