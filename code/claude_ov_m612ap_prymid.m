//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_prymid.m     lane 9 ([6,12])   2026-07-25
//
// The whole Abel-Prym sieve rests on ONE inherited hypothesis:
//     Pr := (1-iota) J(E8)  is isogenous over Q to  J(D),
//     D : y^2 = -3x^6 + 24x^3 - 75,   rank J(D)(Q) = 1, J(D)(Q)_tors = 0
// (code/contact6_m612_prym_rank_verifier.m).  That identification came out
// of the bigonal construction and was never certified independently.
//
// This script tests it the cheap decisive way: for every good prime q it
// compares the FULL L-polynomials,
//     L_{E8}(T)  =?  L_{E4}(T) * L_D(T),
// which is exactly "Frobenius on Pr = Frobenius on J(D)" at q.  Agreement
// at many q is the standard practical isogeny certificate (Faltings/Tate:
// equality of all L_q forces isogeny).
//
// Usage: code/claude_magma_slot.sh -b qmax:=101 code/claude_ov_m612ap_prymid.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned qmax  then qmax := 101; elif Type(qmax) eq MonStgElt then qmax := StringToInteger(qmax); end if;
if not assigned memgb then memgb := 8;  elif Type(memgb) eq MonStgElt then memgb := StringToInteger(memgb); end if;
if not assigned nfork then nfork := 8;  elif Type(nfork) eq MonStgElt then nfork := StringToInteger(nfork); end if;
if not assigned wdir  then wdir := "results/prymid"; end if;
SetMemoryLimit(memgb*10^9);

QQ := Rationals(); Px<x> := PolynomialRing(QQ);
fD := -3*x^6 + 24*x^3 - 75;
printf "Prym candidate D : y^2 = %o\n", fD;

writelines := procedure(fn, L)
  fh := Open(fn, "w");
  for s in L do Puts(fh, s); end for;
  Flush(fh); delete fh;
end procedure;

worker := function(q)
  L := [];
  kq := GF(q);
  Fx<X> := RationalFunctionField(kq);
  Py<Y> := PolynomialRing(Fx);
  Q8 := Y^8 + (216*X^4+72*X^3-24*X^2)*Y^4
      + (-1296*X^6-1728*X^5-432*X^4+64*X^3)*Y^2
      + (-3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4);
  Q4 := Y^4 + (216*X^2+72*X-24)*Y^2 + (-1296*X^3-1728*X^2-432*X+64)*Y
      + (-3888*X^4-2592*X^3+432*X^2+288*X-48);
  if not IsIrreducible(Q8) then return ["SKIP reducible E8"]; end if;
  F8 := FunctionField(Q8);
  if Genus(F8) ne 4 then return [Sprintf("SKIP E8 genus %o", Genus(F8))]; end if;
  F4 := FunctionField(Q4);
  if Genus(F4) ne 2 then return [Sprintf("SKIP E4 genus %o", Genus(F4))]; end if;
  fq := PolynomialRing(kq) ! [kq!c : c in Coefficients(fD)];
  if Discriminant(fq) eq 0 then return ["SKIP D bad reduction"]; end if;
  CD := HyperellipticCurve(fq);
  if Genus(CD) ne 2 then return ["SKIP D genus"]; end if;
  L8 := LPolynomial(F8);
  L4 := LPolynomial(F4);
  LD := LPolynomial(CD);
  Append(~L, Sprintf("L8 %o", Coefficients(L8)));
  Append(~L, Sprintf("L4 %o", Coefficients(L4)));
  Append(~L, Sprintf("LD %o", Coefficients(LD)));
  ok := (L8 eq L4*LD);
  Append(~L, Sprintf("MATCH %o", ok));
  if not ok then
    // also try the quadratic twist / the other sign conventions
    Append(~L, Sprintf("QUOT %o", L8 div L4 eq LD));
  end if;
  Append(~L, Sprintf("ORD8 %o", Evaluate(L8,1)));
  Append(~L, Sprintf("ORD4 %o", Evaluate(L4,1)));
  Append(~L, Sprintf("ORDD %o", Evaluate(LD,1)));
  return L;
end function;

if assigned single then
  s := Type(single) eq MonStgElt select StringToInteger(single) else single;
  for t in worker(s) do printf "q=%o %o\n", s, t; end for;
  print "PRYMID_SINGLE_DONE"; quit;
end if;

primes := [q : q in PrimesInInterval(7, qmax)];
printf "primes: %o\n", primes;
running := 0;
for q in primes do
  if Fork() eq 0 then
    tq := Cputime();
    L := [];
    try L := worker(q); catch err L := ["SKIP error " cat Sprint(err`Object)]; end try;
    Append(~L, Sprintf("TIME %o", Cputime(tq)));
    writelines(wdir cat Sprintf("/prymid_q%o.txt", q), L);
    quit;
  end if;
  running +:= 1;
  if running mod nfork eq 0 then WaitForAllChildren(); end if;
end for;
WaitForAllChildren();
print "PRYMID_DONE";
quit;
