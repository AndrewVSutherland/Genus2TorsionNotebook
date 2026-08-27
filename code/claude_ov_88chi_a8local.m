// claude_ov_88chi_a8local.m -- Lane 5 (2026-07-25).  Extended chart-level local
// viability probe for [8,8] on the A(8) = M_1(8,4) chart (r,p,t).
//
// Positive control: reproduce the 2026-07-18 counts contains88 = 6/98 (q=7),
// 39/589 (q=11), 106/1097 (q=13)  [code/fable_a12_a8_localprobe.m].
// Then extend to q = 17 .. 97 to see whether the local picture stays healthy
// (a degradation would be the first evidence of a real obstruction).
//
// Params: QQ (the prime), PART, NPARTS  (shard over the first coordinate r).
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned QQ then printf "need QQ\n"; quit; end if;
if Type(QQ) eq MonStgElt then QQ := StringToInteger(QQ); end if;
if not assigned PART then PART := 0; elif Type(PART) eq MonStgElt then PART := StringToInteger(PART); end if;
if not assigned NPARTS then NPARTS := 1; elif Type(NPARTS) eq MonStgElt then NPARTS := StringToInteger(NPARTS); end if;

function ContainsInv(iv, ta)
    k := #iv; m := #ta;
    if k lt m then return false; end if;
    return &and[ IsDivisibleBy(iv[k-m+j], ta[j]) : j in [1..m] ];
end function;

q := QQ;
K := GF(q);
Pol<x> := PolynomialRing(K);
els := [ K | i : i in [0..q-1] ];
tested := 0; good := 0; c8 := 0; c44 := 0; c84 := 0; c88 := 0; c168 := 0;
ex := [];
idx := 0;
for rv in els do
  idx +:= 1;
  if (idx mod NPARTS) ne PART then continue; end if;
  if rv eq 0 then continue; end if;
  for pv in els do
    for tv in els do
      if tv eq 0 then continue; end if;
      e := tv^2 - 2*pv*tv/rv;
      d := e + 2*pv - rv^2;
      lam := rv/tv;
      aa := rv^2 - lam;
      bb := 2*rv*pv - 2*lam*(pv + rv*tv) + 2*rv*lam;
      cc := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
            - lam*(rv^2 + e) + 2*lam*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
      tested +:= 1;
      Qc := x^2 + d;
      qc := aa*x^2 + bb*x + cc;
      f := qc*(Qc^2 + qc);
      if Degree(f) ne 6 or not IsSeparable(f) then continue; end if;
      okc := true;
      try
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        n := #J;
        if n mod 64 ne 0 then
          good +:= 1;
          G := AbelianGroup(J);
          iv := Invariants(G);
          if ContainsInv(iv,[8]) then c8 +:= 1; end if;
          if ContainsInv(iv,[4,4]) then c44 +:= 1; end if;
          continue;
        end if;
        G := AbelianGroup(J);
      catch ee okc := false; end try;
      if not okc then continue; end if;
      good +:= 1;
      iv := Invariants(G);
      if ContainsInv(iv, [8])   then c8   +:= 1; end if;
      if ContainsInv(iv, [4,4]) then c44  +:= 1; end if;
      if ContainsInv(iv, [4,8]) then c84  +:= 1; end if;
      if ContainsInv(iv, [8,8]) then
        c88 +:= 1;
        if #ex lt 3 then Append(~ex, <rv,pv,tv,iv>); end if;
      end if;
      if ContainsInv(iv, [8,16]) then c168 +:= 1; end if;
    end for;
  end for;
end for;
printf "A8LOCAL q=%o part=%o tested=%o good=%o c8=%o c44=%o c48=%o c88=%o c816=%o\n",
       q, PART, tested, good, c8, c44, c84, c88, c168;
for e in ex do printf "A8HIT q=%o r=%o p=%o t=%o inv=%o\n", q, e[1], e[2], e[3], e[4]; end for;
printf "A8LOCAL_DONE q=%o part=%o\n", q, PART;
quit;
