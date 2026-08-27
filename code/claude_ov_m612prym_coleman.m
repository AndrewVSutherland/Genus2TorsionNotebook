//////////////////////////////////////////////////////////////////////
// claude_ov_m612prym_coleman.m   (lane 9, overnight 2026-07-25)
//
// PRECISION-SAFE, FORKED Prym-Chabauty driver for E8, the [6,12]
// gatekeeper.  Replaces the WITHDRAWN
// code/contact6_m612_E8_coleman_server.m, whose defects were:
//   (B1) annihilator coefficients known only to O(p^4) were coerced into
//        a precision-N parent, so Tuitman's zeros_on_disk (which reads
//        Nv = Precision(Parent(v[1][1]))) rounded N-4 fabricated digits
//        into the root-finding polynomial;
//   (B2) the normalisation hard-coded  37^mv  instead of  p^mv  (so any
//        run at p <> 37 -- e.g. the README's suggested p := 79 -- would
//        silently produce garbage);
//   (B3) the post-processing DISCARDED each zero's Tuitman `inf` flag, so
//        a zero in an infinity disk was rationally reconstructed and
//        exact-fibre-tested on the WRONG fibre (x instead of 1/x);
//   (B4) no per-disk precision accounting at all, and no TOTAL line.
//
// This driver:
//   * takes p and N from the command line and uses p^mv everywhere;
//   * reports AbsolutePrecision at every stage and asserts the transport
//     identities to the working precision;
//   * re-implements the per-disk root finding (instead of calling
//     zeros_on_disk) so that the working precision Nv of each disk is
//        Nv = Min( Nann, NII, min_j ( j + AbsolutePrecision(f_j) ) )
//     where f = sum f_j z^j is the Chabauty series on the disk and Nann
//     is the number of reliable digits of the annihilator; every
//     coefficient of the root-finding polynomial h_j = p^j f_j is then
//     known to at least Nv digits, so NO fabricated digit ever enters
//     my_roots_Zpt;
//   * retains P`inf for every zero and inverts the coordinate (e = 1/x)
//     in infinity disks before reconstruction and the exact fibre test;
//   * forks one child per residue disk (Magma 2.29-8 POSIX Fork), doing
//     coleman_data / Q_points / the transport ONCE in the parent.
//
// Tuitman library pinned at commit 9e49ae62e3e3bd7e3a6591266595be965cec1a97
// (github.com/jtuitman/Coleman).  Run from inside the Coleman directory.
//
// Command line:
//   magma -b p:=37 N:=28 mode:=pair  dir:=~/run_N28
//   magma -b p:=37 N:=28 mode:=disks dir:=~/run_N28 \
//         a3:=... n3:=... a4:=... n4:=...
//   magma -b p:=37 N:=28 mode:=both  dir:=~/run_N28
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "coleman.m";

Q8 := y^8 + (216*x^4+72*x^3-24*x^2)*y^4 + (-1296*x^6-1728*x^5-432*x^4+64*x^3)*y^2
     + (-3888*x^8-2592*x^7+432*x^6+288*x^5-48*x^4);

// ------------------------------- parameters
if not assigned p     then p := 37;     else p := StringToInteger(p);         end if;
if not assigned N     then N := 20;     else N := StringToInteger(N);         end if;
if not assigned epts  then epts := 100; else epts := StringToInteger(epts);   end if;
if not assigned ee    then ee := 300;   else ee := StringToInteger(ee);       end if;
// eeb: the expansion length for the boundary-boundary integral only
if not assigned eeb   then eeb := ee;   else eeb := StringToInteger(eeb);      end if;
// skip7 = 1 uses the THEOREM Int_{b1}^{b2} omega = 0 (both boundary places are
// iota-fixed and omega is anti-invariant) instead of computing that integral,
// which is the single most expensive and least precise one.
if not assigned skip7 then skip7 := 0;  else skip7 := StringToInteger(skip7);  end if;
if not assigned mode  then mode := "both"; end if;
if not assigned dir   then dir := "/tmp/m612prym"; end if;
// NOTE: the command-line names are deliberately prefixed (pa3, ...) because
// loading Tuitman's coleman.m leaves short identifiers such as a3 assigned
// at top level, which makes a bare "assigned a3" test unreliable.
havepair := false;
if assigned pa3 and assigned pn3 and assigned pa4 and assigned pn4 then
  a3 := StringToInteger(pa3); n3 := StringToInteger(pn3);
  a4 := StringToInteger(pa4); n4 := StringToInteger(pn4);
  havepair := true;
end if;
System("mkdir -p " cat dir);

printf "=== claude_ov_m612prym_coleman ===\n";
printf "p = %o, N = %o, epts = %o, ee = %o, mode = %o, dir = %o\n",
    p, N, epts, ee, mode, dir;

tstart := Cputime();
data := coleman_data(Q8, p, N);
printf "coleman_data ready: g = %o, N = %o  (%o s)\n", data`g, N, Cputime(tstart);

Qp := pAdicField(p, N);
mkpad := function(a, n)          // a mod p^n as an element of Qp with abs prec n
  return (Qp!(a mod p^n)) + O(Qp!p^n);
end function;
padout := function(u)            // "lift absprec" for a p-adic of valuation >= 0
  n := AbsolutePrecision(u);
  assert Valuation(u) ge 0;
  return IntegerToString(Integers()!(pAdicRing(p,n)!u)) cat " " cat IntegerToString(n);
end function;

evalQ8 := func<X0, Y0 | &+[ Evaluate(ChangeRing(Coefficient(Q8, i), Parent(X0)), X0)*Y0^i
                            : i in [0..8]]>;

qpts := Q_points(data, 1000);
printf "rational points found by Q_points: %o\n", #qpts;
for i in [1..#qpts] do
  printf "  qpt %o: x = %o, inf = %o\n", i, qpts[i]`x, qpts[i]`inf;
end for;

// ------------------------------- stage 2 : the pairing
// Generator of J(C2'_min)(Q) (rank 1): Mumford (x^2+2x+4, 5x+5);
// support x = -1 +- sqrt(-3), y = +-5 sqrt(-3).  To C2'par: t = 2x+1,
// s = y/3.  Bigonal transport: xE4 = -72/(t^2+3), v = -72t/(t^2+3),
// u = 6(xE4^2-6xE4-36)v/(xE4^2(xE4+24)); w,w' roots of T^2-sT+u;
// Tuitman coords X = e = (xE4/3+2)/(w^2-xE4), Y = e*w.

if not havepair then
  Pq<T> := PolynomialRing(Qp);
  rts := [r[1] : r in Roots(T^2+T+1)];
  resid := [Integers()!(GF(p)!(Integers()!(pAdicRing(p,1)!r))) : r in rts];
  if resid[1] le resid[2] then z3 := rts[1]; else z3 := rts[2]; end if;
  lam := 2*z3+1;          // sqrt(-3)
  printf "zeta_3 residues %o (chose %o), v(lam^2+3) = %o\n",
      resid, Minimum(resid), Valuation(lam^2+3);
  assert Valuation(lam^2 + 3) ge N-2;

  xm := [-1+lam, -1-lam];
  ts := [2*xm[1]+1, 2*xm[2]+1];
  ss := [(5*xm[1]+5)/3, (5*xm[2]+5)/3];

  Kt<tt> := FunctionField(RationalField());
  xtf := -72/(tt^2+3); vtf := -72*tt/(tt^2+3);
  utf := 6*(xtf^2-6*xtf-36)*vtf/(xtf^2*(xtf+24));
  Atf := (xtf^4/4 - 15*xtf^3/2 + 324*xtf + 972)/(xtf^3*(xtf+24)/24);
  Stf := 2*Atf + 2*utf;
  numS := Numerator(Stf); denS := Denominator(Stf);
  f2t := numS*denS;
  f2sft := &*[ fe[1]^(fe[2] mod 2) : fe in Factorization(f2t) ] * LeadingCoefficient(f2t);
  sq2 := f2t / f2sft;
  assert Denominator(sq2) eq 1;
  sqpoly := Sqrt(Numerator(sq2));
  printf "den deg %o, sq deg %o, f2sf deg %o\n",
      Degree(denS), Degree(sqpoly), Degree(Numerator(f2sft));

  pts_D1 := [];
  for i in [1,2] do
    t0 := ts[i]; s0m := ss[i];
    sqv  := &+[Qp | Qp!Coefficient(sqpoly,k)*t0^k : k in [0..Degree(sqpoly)]];
    denv := &+[Qp | Qp!Coefficient(denS,k)*t0^k : k in [0..Degree(denS)]];
    s0 := s0m * sqv/denv;
    x0 := -72/(t0^2+3);
    v0 := -72*t0/(t0^2+3);
    u0 := 6*(x0^2-6*x0-36)*v0/(x0^2*(x0+24));
    A0 := (x0^4/4 - 15*x0^3/2 + 324*x0 + 972)/(x0^3*(x0+24)/24);
    cons := s0^2 - 2*A0 - 2*u0;
    printf "  transport %o: v(s^2-(2A+2u)) = %o (abs prec %o)\n",
        i, Valuation(cons), AbsolutePrecision(cons);
    assert Valuation(cons) ge AbsolutePrecision(cons) - 2;
    dsc := s0^2 - 4*u0;
    sq := Sqrt(dsc);
    for sgn in [1,-1] do
      w0 := (s0 + sgn*sq)/2;
      e0 := (x0/3 + 2)/(w0^2 - x0);
      X0 := e0; Y0 := e0*w0;
      err := evalQ8(X0, Y0);
      printf "  D1 point %o/%o: v(X) = %o, v(model err) = %o (abs prec %o)\n",
          i, sgn, Valuation(X0), Valuation(err), AbsolutePrecision(err);
      assert Valuation(err) ge AbsolutePrecision(err) - 2;
      Append(~pts_D1, <X0, Y0>);
    end for;
  end for;

  pts_Dinf_fin := [ <Qp!(-2)/3, 2*lam/3>, <Qp!(-2)/3, -2*lam/3> ];
  for P in pts_Dinf_fin do
    err := evalQ8(P[1], P[2]);
    printf "  Dinf finite pt: v(model err) = %o (abs prec %o)\n",
        Valuation(err), AbsolutePrecision(err);
    assert Valuation(err) ge AbsolutePrecision(err) - 2;
  end for;

  bb := qpts[1];
  rpoly := data`r;
  W0m := data`W0;
  allpts := [];
  for P in pts_D1 cat pts_Dinf_fin do
    X0 := P[1]; Y0 := P[2];
    rv := &+[Qp | Qp!Coefficient(rpoly,k)*X0^k : k in [0..Degree(rpoly)]];
    if Valuation(rv) eq 0 then
      Append(~allpts, set_point(X0, Y0, data));
    else
      bvec := [];
      for j in [1..8] do
        bj := Qp!0;
        for k in [1..8] do
          w0jk := W0m[j,k];
          if w0jk ne 0 then
            nn := Numerator(w0jk); dd := Denominator(w0jk);
            nv := &+[Qp | Qp!Coefficient(nn,h)*X0^h : h in [0..Degree(nn)]];
            dv := &+[Qp | Qp!Coefficient(dd,h)*X0^h : h in [0..Degree(dd)]];
            bj +:= (nv/dv) * Y0^(k-1);
          end if;
        end for;
        Append(~bvec, bj);
      end for;
      printf "  bad-disk point (v(x) = %o): set_bad_point\n", Valuation(X0);
      Append(~allpts, set_bad_point(X0, bvec, false, data));
    end if;
  end for;
  targets := allpts cat [qpts[2]];
  eeach   := [epts, epts, epts, epts, epts, epts, eeb];
  if skip7 eq 1 then
    targets := allpts;  eeach := [epts, epts, epts, epts, epts, epts];
    printf "skip7: using the theorem Int_{b1}^{b2} omega_anti = 0 (both boundary places are iota-fixed)\n";
  end if;

  // ---- fork one child per Coleman integral (7 of them)
  for c in [1..#targets] do
    pid := Fork();
    if pid eq 0 then
      SetMemoryLimit(12*10^9);
      SetOutputFile(dir cat "/stdout_pairint" cat IntegerToString(c) cat ".txt" : Overwrite := true);
      F := Open(dir cat "/pairint" cat IntegerToString(c) cat ".txt", "w");
      ti := Cputime();
      II, NII := coleman_integrals_on_basis(bb, targets[c], data : e := eeach[c]);
      Puts(F, padout(II[3]));
      Puts(F, padout(II[4]));
      Puts(F, IntegerToString(NII) cat " " cat IntegerToString(Round(Cputime(ti))));
      Flush(F); delete F;
      quit;
    end if;
  end for;
  WaitForAllChildren();

  Ints3 := []; Ints4 := []; NIIs := [];
  for c in [1..#targets] do
    L := Split(Read(dir cat "/pairint" cat IntegerToString(c) cat ".txt"), "\n");
    s1 := Split(L[1], " "); s2 := Split(L[2], " "); s3 := Split(L[3], " ");
    Append(~Ints3, mkpad(StringToInteger(s1[1]), StringToInteger(s1[2])));
    Append(~Ints4, mkpad(StringToInteger(s2[1]), StringToInteger(s2[2])));
    Append(~NIIs, StringToInteger(s3[1]));
    printf "  integral %o: I3 = %o, I4 = %o ; NII = %o (%o s)\n",
        c, Ints3[c], Ints4[c], NIIs[c], s3[2];
  end for;

  // ---- validations that do not depend on the transport being right
  if skip7 eq 1 then
    Append(~Ints3, Qp!0); Append(~Ints4, Qp!0);
  else
    printf "  VALIDATION boundary anti-invariance (b -> b2, must vanish): v(I3) = %o (abs prec %o), v(I4) = %o (abs prec %o)\n",
        Valuation(Ints3[7]), AbsolutePrecision(Ints3[7]),
        Valuation(Ints4[7]), AbsolutePrecision(Ints4[7]);
  end if;
  printf "  VALIDATION Dinf-pair cancellation: v(I5+I6)_3 = %o , v(I5+I6)_4 = %o\n",
      Valuation(Ints3[5]+Ints3[6]), Valuation(Ints4[5]+Ints4[6]);

  pair3 := &+[Ints3[i] : i in [1..4]] - Ints3[5] - Ints3[6] - Ints3[7];
  pair4 := &+[Ints4[i] : i in [1..4]] - Ints4[5] - Ints4[6] - Ints4[7];
  printf "PAIRING basis[3] = %o  (val %o, abs prec %o)\n",
      pair3, Valuation(pair3), AbsolutePrecision(pair3);
  printf "PAIRING basis[4] = %o  (val %o, abs prec %o)\n",
      pair4, Valuation(pair4), AbsolutePrecision(pair4);
  n3 := AbsolutePrecision(pair3); n4 := AbsolutePrecision(pair4);
  a3 := Integers()!(pAdicRing(p,n3)!pair3);
  a4 := Integers()!(pAdicRing(p,n4)!pair4);
  printf "PAIRDATA a3:=%o n3:=%o a4:=%o n4:=%o\n", a3, n3, a4, n4;
else
  pair3 := mkpad(a3, n3);
  pair4 := mkpad(a4, n4);
  printf "PAIRING supplied: pair3 = %o (abs prec %o), pair4 = %o (abs prec %o)\n",
      pair3, AbsolutePrecision(pair3), pair4, AbsolutePrecision(pair4);
end if;

mv := Min(Valuation(pair3), Valuation(pair4));
printf "mv = %o (normalisation by p^mv with p = %o)\n", mv, p;
c3 :=  pair4/p^mv;
c4 := -pair3/p^mv;
Nann := Min(AbsolutePrecision(c3), AbsolutePrecision(c4));
printf "omega_ann = (%o)*basis[3] + (%o)*basis[4]\n", c3, c4;
printf "ANNIHILATOR RELIABLE DIGITS Nann = %o\n", Nann;
if Nann lt 12 then
  printf "WARNING: fewer than 12 reliable p-adic digits in the annihilator; raise N.\n";
end if;

if mode eq "pair" then
  print "M612PRYM_PAIR_DONE";
  quit;
end if;

// ------------------------------- stage 3 : the disk sweep, one child per disk
for i := 1 to #qpts do
  _, index := local_data(qpts[i], data);
  data := update_minpolys(data, qpts[i]`inf, index);
  if is_bad(qpts[i], data) then
    if is_very_bad(qpts[i], data) then
      xt, bt, index := local_coord(qpts[i], tadicprec(data, ee), data);
      qpts[i]`xt := xt; qpts[i]`bt := bt; qpts[i]`index := index;
    end if;
  else
    xt, bt, index := local_coord(qpts[i], tadicprec(data, 1), data);
    qpts[i]`xt := xt; qpts[i]`bt := bt; qpts[i]`index := index;
  end if;
end for;
Qppoints, data := Qp_points(data : points := qpts);
ND := #Qppoints;
printf "residue disks: %o\n", ND;
bb := qpts[1];

for i := 1 to ND do
  pid := Fork();
  if pid eq 0 then
    SetMemoryLimit(12*10^9);
    SetOutputFile(dir cat "/stdout_disk" cat IntegerToString(i) cat ".txt" : Overwrite := true);
    F := Open(dir cat "/disk" cat IntegerToString(i) cat ".txt", "w");
    tdisk := Cputime();
    QQ := Rationals(); Pw<W> := PolynomialRing(QQ);
    P2 := Qppoints[i];
    _, index := local_data(P2, data);
    data := update_minpolys(data, P2`inf, index);
    if is_bad(P2, data) then
      xt, bt, index := local_coord(P2, tadicprec(data, ee), data);
    else
      xt, bt, index := local_coord(P2, tadicprec(data, 1), data);
    end if;
    P2`xt := xt; P2`bt := bt; P2`index := index;
    II, NII := coleman_integrals_on_basis(bb, P2, data : e := ee);
    tinyz, xts, bts, Ntiny := tiny_integrals_on_basis_to_z(P2, data);
    ser := c3*(II[3] + tinyz[3]) + c4*(II[4] + tinyz[4]);
    dg := Degree(ser);
    precs := [ j + AbsolutePrecision(Coefficient(ser,j)) : j in [0..dg] ];
    Nv := Max(Min([Nann, NII] cat precs), 1);
    f0 := Coefficient(ser, 0);
    Puts(F, Sprintf("DISKINFO %o inf=%o NII=%o Nann=%o minjprec=%o Nv=%o vf0=%o deg=%o",
        i, P2`inf, NII, Nann, Min(precs), Nv, Valuation(f0), dg));
    Zp := pAdicRing(p, Nv);
    Zpt<Z> := PolynomialRing(Zp);
    h := Zpt!0;
    for j := 0 to dg do
      cf := Coefficient(ser, j);
      if cf eq 0 then continue; end if;
      h +:= (Integers()!(p^j*(RationalField()!cf)))*Z^j;
    end for;
    if h eq 0 then
      Puts(F, Sprintf("DISK %o ZEROFUNCTION", i));
      Puts(F, Sprintf("DISKDONE %o -1 %o", i, Round(Cputime(tdisk))));
      Flush(F); delete F; quit;
    end if;
    Puts(F, Sprintf("NEWTON %o %o", i,
        [ <j, Valuation(Coefficient(h,j))> : j in [0..Min(dg,4)] ]));
    zeros := my_roots_Zpt(h);
    Puts(F, Sprintf("DISKDONE %o %o %o", i, #zeros, Round(Cputime(tdisk))));
    for zz in zeros do
      z := zz[1]; Nz := zz[2];
      x0 := Evaluate(xts, p*z);
      bvec := Eltseq(Evaluate(bts, p*z));
      P := set_bad_point(x0, bvec, P2`inf, data);
      isinf := P`inf;
      xv := P`x;
      Puts(F, Sprintf("ZERO %o x=%o inf=%o absprec=%o rootprec=%o",
          i, xv, isinf, AbsolutePrecision(xv), Nz));
      skip := false;
      if isinf then
        if Valuation(xv) ge AbsolutePrecision(xv) then
          Puts(F, Sprintf("CLASS %o INFINITY", i)); skip := true;
        else
          e0p := 1/xv;
        end if;
      else
        e0p := xv;
      end if;
      if not skip then
        vv := Minimum(Valuation(e0p), 0);
        kk := Minimum(AbsolutePrecision(e0p) - vv, Nv);
        if kk le 1 then
          Puts(F, Sprintf("CLASS %o NOPRECISION", i));
        else
          num := e0p*p^(-vv);
          ii := Integers()!(pAdicRing(p,kk)!num);
          ok, q0 := RationalReconstruction(Integers(p^kk)!ii);
          if ok then q0 := q0*p^vv; end if;
          if not ok then
            Puts(F, Sprintf("CLASS %o MOCK_NORECON kk=%o", i, kk));
          elif q0 eq 0 then
            Puts(F, Sprintf("CLASS %o BOUNDARY e=0", i));
          else
            fib := Pw ! [Evaluate(Coefficient(Q8, j), q0) : j in [0..8]];
            rr := Roots(fib);
            if #rr gt 0 then
              Puts(F, Sprintf("CLASS %o RATIONALPOINT e=%o y=%o", i, q0, [r[1] : r in rr]));
            else
              Puts(F, Sprintf("CLASS %o MOCK_EMPTYFIBRE e=%o kk=%o", i, q0, kk));
            end if;
          end if;
        end if;
      end if;
    end for;
    Flush(F); delete F;
    quit;
  end if;
end for;
WaitForAllChildren();

totzero := 0; nfail := 0;
for i := 1 to ND do
  fn := dir cat "/disk" cat IntegerToString(i) cat ".txt";
  ok := true;
  try
    txt := Read(fn);
  catch err
    ok := false;
  end try;
  if not ok then
    printf "DISK %o : NO OUTPUT (child failed)\n", i; nfail +:= 1; continue;
  end if;
  for l in Split(txt, "\n") do
    if #l gt 0 then printf "%o\n", l; end if;
    if #l ge 8 and l[1..8] eq "DISKDONE" then
      s := Split(l, " ");
      nz := StringToInteger(s[3]);
      if nz ge 0 then totzero +:= nz; else nfail +:= 1; end if;
    end if;
  end for;
end for;
printf "TOTAL ZEROS over %o residue disks = %o  (failed disks: %o)\n", ND, totzero, nfail;
print "M612PRYM_COLEMAN_DONE";
quit;
