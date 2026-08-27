// claude_ov_88verify.m -- Lane 5 (2026-07-25, third session).
//
// UNCONDITIONAL re-verification of the [8,8] negative on the Nicholls Lambda_334
// double-stage-1 (DS1) locus.
//
// Why this script exists.  The 108-base conic harvest
// (results/claude_ov_88conic2_partial_harvest.log) was produced by a Magma session
// that had SetClassGroupBounds("GRH") set globally, because the elliptic
// MordellWeilGroup route needed it.  Everything downstream of the *member*
// construction (exact TorsionSubgroup, the x-T gate) is class-group free, but that
// is an argument, not a measurement.  This script re-derives, in a fresh session
// with NO GRH bound set, for every DISTINCT curve of the harvest:
//
//   * the member (s,t,v) -> integral sextic  fint     (deterministic, no search)
//   * the irreducible factor shape of fint
//   * the EXACT rational torsion  J(Q)_tors = TorsionSubgroup(J)
//   * for each order-4 class D:  GROUND TRUTH  "D in 2*J(Q)_tors"  (exact -- any
//     rational half of a torsion point is torsion, so 2*T is the whole of
//     2J(Q) cap T), and the x-T GATE (common-lambda test, the same code as
//     claude_ov_88gate_validate.m), and the two square-class defects
//     def_i = lambda_ref*lambda_i that the conic layer reports.
//
// [8,8] requires TWO cyclic factors divisible by 8; on this locus the 2-torsion
// has rank 2, so J(Q)_tors contains Z/8 x Z/8 <=> ALL order-4 classes are halvable.
//
// POSITIVE CONTROL: with CONTROL:=1 the same Gate/ground-truth code is run on the
// LMFDB curves in data/claude_ov_88gate_curves.txt, so the log contains, in the
// same format, curves on which the gate DOES return a nonzero count.
//
// Params: FILE (default data/claude_ov_88verify_stv.txt, lines "s t v tag"),
//         CONTROL (default 1), CFILE (default data/claude_ov_88gate_curves.txt),
//         NCH (children, default 8), OUTDIR
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);  Z := Integers();

if not assigned FILE then FILE := "data/claude_ov_88verify_stv.txt"; end if;
if not assigned CFILE then CFILE := "data/claude_ov_88gate_curves.txt"; end if;
if not assigned CONTROL then CONTROL := 1; elif Type(CONTROL) eq MonStgElt then CONTROL := StringToInteger(CONTROL); end if;
if not assigned NCH then NCH := 8; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned OUTDIR then
  OUTDIR := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/verify";
end if;
System("mkdir -p " cat OUTDIR);

function SqFree(r)
  if r eq 0 then return 0; end if;
  m := Numerator(r)*Denominator(r);
  s := Sign(m);  m := AbsoluteValue(m);
  res := s;
  for pf in Factorization(m) do if IsOdd(pf[2]) then res *:= pf[1]; end if; end for;
  return res;
end function;

// member of Lambda_334 from (s,t,v); identical formulas to claude_ov_88conic_base2.m
function MemberSTV(sv, tv, vv)
  if sv eq 0 or tv eq 0 or tv^2 eq 1 then return false, _, _; end if;
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 then return false, _, _; end if;
  den := -sv^2*tv*Av*vv^2 + tv;
  if den eq 0 then return false, _, _; end if;
  uv := (-sv^2*Av*vv^2 - 2*Av*vv - 1)/den;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false, _, _; end if;
  av := Av/(1 - tv^2);
  bv := Av/(uv^2*sv^2 + 1 - tv^2);
  cv := tv^2;
  d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)
          *(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
  if d2 eq 0 then return false, _, _; end if;
  l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
  l2 := -x^2 + bv*cv;
  l3 := x^2 - av;
  f1 := l1*l2*l3;
  if Degree(f1) ne 6 then return false, _, _; end if;
  f := d2*f1/LeadingCoefficient(f1);
  if Discriminant(f) eq 0 then return false, _, _; end if;
  dd := LCM([Denominator(cc) : cc in Coefficients(f)]);
  fint := P![ cc*dd^2 : cc in Coefficients(f) ];
  cont := GCD([Z!cc : cc in Coefficients(fint)]);
  sqf := 1;
  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  fint := P![ (Z!cc) div sqf : cc in Coefficients(fint) ];
  return true, fint, [av,bv,cv,Av,uv];
end function;

// The GATE: does there exist lambda in Q^* with u_D in lambda*(K_i^*)^2 for all i?
// returns status, lams, ds, compat
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
      if not isq then return "normfail", [], [], false; end if;
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

// full analysis of one integral sextic/quintic; returns a list of output strings
function Analyse(tag, fint)
  out := [];
  C := HyperellipticCurve(fint);  J := Jacobian(C);
  T, mT := TorsionSubgroup(J);
  inv := Invariants(T);
  facs := [ g[1]/LeadingCoefficient(g[1]) : g in Factorization(fint) ];
  shape := Sort([ Degree(g) : g in facs ]);
  dbl := { mT(2*g) : g in T };
  n4 := 0; ntruth := 0; ngate := 0; nnorm := 0; ndegu := 0; ndis := 0;
  lines := [];
  for g in T do
    if Order(g) ne 4 then continue; end if;
    n4 +:= 1;
    D := mT(g);  uD := D[1];
    tr := D in dbl;
    if tr then ntruth +:= 1; end if;
    if Degree(uD) ne D[3] then ndegu +:= 1; end if;
    st, lams, ds, cp := Gate(uD, facs);
    if st eq "normfail" then nnorm +:= 1; end if;
    if cp then ngate +:= 1; end if;
    if cp ne tr then
      ndis +:= 1;
      Append(~lines, Sprintf("  DISAGREE %o u=%o st=%o gate=%o truth=%o lam=%o ds=%o", tag, uD, st, cp, tr, lams, ds));
    end if;
    if st eq "ok" then
      // square-class defects relative to the first component
      def := [ SqFree(lams[1]*lams[i]) : i in [2..#lams] ];
      Append(~lines, Sprintf("  D4 %o ds=%o lam=%o def=%o gate=%o truth=%o", tag, ds, lams, def, cp, tr));
    else
      Append(~lines, Sprintf("  D4 %o st=%o", tag, st));
    end if;
  end for;
  Append(~out, Sprintf("VERIFY %o deg=%o shape=%o tors=%o order4=%o normfail=%o degu=%o gate_pass=%o truth_pass=%o disagree=%o",
         tag, Degree(fint), shape, inv, n4, nnorm, ndegu, ngate, ntruth, ndis));
  for l in lines do Append(~out, l); end for;
  return out;
end function;

// ---------- job list ----------
jobs := [];   // <kind, tag, data...>
for ln in Split(Read(FILE), "\n") do
  if ln eq "" or ln[1] eq "#" then continue; end if;
  tk := Split(ln, " ");
  if #tk lt 3 then continue; end if;
  Append(~jobs, <"stv", tk[1] cat "|" cat tk[2] cat "|" cat tk[3], tk[1], tk[2], tk[3]>);
end for;
nfam := #jobs;
if CONTROL eq 1 then
  for ln in Split(Read(CFILE), "\n") do
    if ln eq "" or ln[1] eq "#" then continue; end if;
    pr := Split(ln, "|");
    if #pr lt 3 then continue; end if;
    lab := pr[1];  while lab[#lab] eq " " do lab := Substring(lab,1,#lab-1); end while;
    Append(~jobs, <"lmfdb", "CTRL:" cat lab, pr[2], pr[3], "">);
  end for;
end if;
printf "VERIFY_LOADED family=%o control=%o total=%o children=%o\n", nfam, #jobs-nfam, #jobs, NCH;

for c in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then
    F := Open(OUTDIR cat "/part" cat IntegerToString(c) cat ".txt", "w");
    i := c+1;
    while i le #jobs do
      jb := jobs[i];
      try
        if jb[1] eq "stv" then
          ok, fint := MemberSTV(Q!eval(jb[3]), Q!eval(jb[4]), Q!eval(jb[5]));
          if not ok then
            Puts(F, Sprintf("VERIFY %o DEGENERATE", jb[2]));
          else
            for l in Analyse(jb[2], fint) do Puts(F, l); end for;
          end if;
        else
          fc := [ Q!StringToInteger(z) : z in Split(jb[3], ",") ];
          hc := [ Q!StringToInteger(z) : z in Split(jb[4], ",") ];
          FF := 4*P!fc + (P!hc)^2;
          dd := LCM([Denominator(cc) : cc in Coefficients(FF)]);
          FF := P![ Z!(cc*dd^2) : cc in Coefficients(FF) ];
          cont := GCD([Z!cc : cc in Coefficients(FF)]);
          sqf := 1;  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
          FF := P![ (Z!cc) div sqf : cc in Coefficients(FF) ];
          for l in Analyse(jb[2], FF) do Puts(F, l); end for;
        end if;
      catch e
        Puts(F, Sprintf("VERIFY %o ERROR %o", jb[2], e`Object));
      end try;
      Flush(F);
      i +:= NCH;
    end while;
    Puts(F, Sprintf("CHILD_DONE %o", c));
    Flush(F); delete F;
    quit;
  end if;
end for;
WaitForAllChildren();
nl := 0;
for c in [0..NCH-1] do
  for l in Split(Read(OUTDIR cat "/part" cat IntegerToString(c) cat ".txt"), "\n") do
    if #l gt 0 then printf "%o\n", l; nl +:= 1; end if;
  end for;
end for;
printf "VERIFY_DONE jobs=%o lines=%o\n", #jobs, nl;
quit;
