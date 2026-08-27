// claude_ov_88chart_verify.m -- Lane 5 (2026-07-25).
// VERIFY a new, much cleaner rational chart for Nicholls' Lambda_334 family.
//
// Derivation (this session).  In Nicholls Prop 5.9.6 only u^2 enters (a,b,c,d2), and
//   b = A/(s^2 u^2 + 1 - t^2),  A = s^2 - t^4 + t^2,  a = A/(1-t^2),  c = t^2.
// Solving the v-conic for u shows the discriminant is (A t)^2 / b, so
//   "v rational"  <=>  b is a SQUARE.                 (measured: 12825/12825 members)
// Put b = beta^2 and eta := s*beta*u.  Then
//   eta^2 = A - beta^2 (1-t^2) = s^2 + (1-t^2)(t^2 - beta^2).            (*)
// So the family is exactly the affine 3-fold (*) in coordinates (s, t, beta, eta),
// with u = eta/(s*beta), and NO v is needed.
//
// x-T conditions of the three prescribed classes on the quintic
// J2 = Jac(y^2 = x(x-1)(x-a)(x-b)(x-c)) become, using
//   a - b = eta^2/(1-t^2),   a - c = s^2/(1-t^2),   b = beta^2, c = t^2 :
//   class {0,a} in 2J2(Q)  <=>  t^2-1 = square  AND  (t^2-1)(s^2-(t^2-1)^2) = square
//   class {c,oo} in 2J2(Q) <=>  t^2-1 = square  AND  t^2 - beta^2 = square
// The second is the SIMPLE form of the old elliptic condition Psi = v(s^2v+1)(Av+1).
//
// Production chart:  t = (m^2+1)/(2m), alpha = (m^2-1)/(2m)      [t^2-1 = alpha^2]
//                    s = (n^2+alpha^4)/(2n)                      [s^2-alpha^4 = square]
//                    beta = t(1-p^2)/(1+p^2), gamma = 2tp/(1+p^2) [t^2-beta^2 = gamma^2]
//                    eta = w/(1+p^2),  w^2 = s^2(1+p^2)^2 - 4 alpha^2 t^2 p^2   (**)
// so the whole DOUBLE-stage-1 locus is the single quartic condition (**) -- a rational
// point search on an elliptic quartic with the visible point (p,w) = (0,s).
//
// This script: (1) checks the chart reproduces the family, (2) checks the two x-T
// criteria against exact Magma torsion on J2, (3) reproduces the recorded members.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);
Z := Integers();

// old chart: from (s,t,v)
function OldMember(sv, tv, vv)
  if sv eq 0 or tv eq 0 or tv^2 eq 1 then return false, _; end if;
  Av := sv^2 - tv^4 + tv^2;  if Av eq 0 then return false, _; end if;
  den := -sv^2*tv*Av*vv^2 + tv;  if den eq 0 then return false, _; end if;
  uv := (-sv^2*Av*vv^2 - 2*Av*vv - 1)/den;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false, _; end if;
  return true, [sv, tv, uv];
end function;

// new chart: from (s,t,beta,eta) with eta^2 = s^2 + (1-t^2)(t^2-beta^2)
function BuildFromU(sv, tv, uv)
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 or sv eq 0 or tv eq 0 or tv^2 eq 1 then return false,_,_,_,_; end if;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false,_,_,_,_; end if;
  av := Av/(1 - tv^2);
  bv := Av/(uv^2*sv^2 + 1 - tv^2);
  cv := tv^2;
  d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)
          *(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
  if d2 eq 0 then return false,_,_,_,_; end if;
  l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
  l2 := -x^2 + bv*cv;  l3 := x^2 - av;
  f1 := l1*l2*l3;
  if Degree(f1) ne 6 then return false,_,_,_,_; end if;
  f := d2*f1/LeadingCoefficient(f1);
  if Discriminant(f) eq 0 then return false,_,_,_,_; end if;
  dd := LCM([Denominator(cc) : cc in Coefficients(f)]);
  fint := P![ cc*dd^2 : cc in Coefficients(f) ];
  cont := GCD([Z!cc : cc in Coefficients(fint)]);
  sqf := 1; for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  fint := P![ (Z!cc) div sqf : cc in Coefficients(fint) ];
  q5 := x*(x-1)*(x-av)*(x-bv)*(x-cv);
  return true, fint, q5, [av,bv,cv], d2;
end function;

printf "=== TEST 1: chart identity  eta^2 = s^2 + (1-t^2)(t^2-beta^2), u = eta/(s*beta) ===\n";
nb := 0; nok := 0;
for sv in [Q| 2,3,5/2,-3/4,7/3] do for tv in [Q| 5/4,3/2,5/3,2,7/5] do
  for vv in [Q| 1,-1,2,-3,1/2,3/2,-5/2] do
    ok, dat := OldMember(sv,tv,vv);
    if not ok then continue; end if;
    uv := dat[3];
    Av := sv^2-tv^4+tv^2;
    bq := Av/(uv^2*sv^2+1-tv^2);
    nb +:= 1;
    isb, bet := IsSquare(bq);
    if not isb then printf "  FAIL b not square s=%o t=%o v=%o\n", sv,tv,vv; continue; end if;
    eta := sv*bet*uv;
    if eta^2 eq sv^2 + (1-tv^2)*(tv^2 - bet^2) then nok +:= 1;
    else printf "  FAIL identity s=%o t=%o v=%o\n", sv,tv,vv; end if;
  end for;
end for; end for;
printf "TEST1 members=%o identity_ok=%o\n", nb, nok;

printf "=== TEST 2: the two stage-1 x-T criteria vs exact Magma torsion on J2 ===\n";
// criterion:  {0,a} divisible <=> t^2-1 sq and (t^2-1)(s^2-(t^2-1)^2) sq
//             {c,oo} divisible <=> t^2-1 sq and t^2-beta^2 sq
nt := 0; agree0a := 0; agreecoo := 0;
for sv in [Q| 2,3,5/2,-3/4,7/3,265/54,4,1/3] do for tv in [Q| 5/4,3/2,5/3,2,7/5,3] do
  for vv in [Q| 1,-1,2,-3,1/2] do
    ok, dat := OldMember(sv,tv,vv);
    if not ok then continue; end if;
    uv := dat[3];
    ok2, fint, q5, abc, d2 := BuildFromU(sv,tv,uv);
    if not ok2 then continue; end if;
    av := abc[1]; bv := abc[2]; cv := abc[3];
    isb, bet := IsSquare(bv);
    if not isb then continue; end if;
    // predicted
    p0a := IsSquare(tv^2-1) and IsSquare((tv^2-1)*(sv^2-(tv^2-1)^2));
    pcoo := IsSquare(tv^2-1) and IsSquare(tv^2-bet^2);
    // actual: build J2 on an INTEGRAL model  Y^2 = d*X(X-d)(X-da)(X-db)(X-dc)
    // (obtained from x = X/d, y = Y/d^3), and test 2-divisibility of the two classes.
    dsc := LCM([Denominator(av), Denominator(bv), Denominator(cv)]);
    G := dsc*x*(x-dsc)*(x-dsc*av)*(x-dsc*bv)*(x-dsc*cv);
    J2 := Jacobian(HyperellipticCurve(G));
    T2, m2 := TorsionSubgroup(J2);
    // classes as points
    D0a := elt<J2 | x*(x-dsc*av), P!0, 2>;
    Dcoo := elt<J2 | x-dsc*cv, P!0, 1>;
    div0a := false; divcoo := false;
    for g in T2 do
      Pg := m2(g);
      if 2*Pg eq D0a then div0a := true; end if;
      if 2*Pg eq Dcoo then divcoo := true; end if;
    end for;
    nt +:= 1;
    if div0a eq p0a then agree0a +:= 1; else
      printf "  MISMATCH {0,a} s=%o t=%o v=%o pred=%o act=%o\n", sv,tv,vv,p0a,div0a; end if;
    if divcoo eq pcoo then agreecoo +:= 1; else
      printf "  MISMATCH {c,oo} s=%o t=%o v=%o pred=%o act=%o beta=%o\n", sv,tv,vv,pcoo,divcoo,bet; end if;
  end for;
end for; end for;
printf "TEST2 members=%o agree_0a=%o agree_coo=%o\n", nt, agree0a, agreecoo;

printf "=== TEST 3: production chart (m,n,p,w) -> member, vs recorded double-stage-1 ===\n";
// recorded: base (m,n) = (3,1/3), v in {-729/17500, 26244/7975, 729/38425}: J2 = [2,2,4,4]
for tup in [ <Q!3, Q!(1/3), Q!(-729/17500)>, <Q!3, Q!(1/3), Q!(26244/7975)>,
             <Q!3, Q!(1/3), Q!(729/38425)>, <Q!2, Q!1, Q!1>, <Q!2, Q!1, Q!(-1)> ] do
  mv := tup[1]; nv := tup[2]; vv := tup[3];
  alv := (mv^2-1)/(2*mv); tv := (mv^2+1)/(2*mv); sv := (nv^2+alv^4)/(2*nv);
  ok, dat := OldMember(sv,tv,vv);
  if not ok then printf "  (m,n,v)=(%o,%o,%o) invalid\n", mv,nv,vv; continue; end if;
  uv := dat[3];
  Av := sv^2-tv^4+tv^2;
  bq := Av/(uv^2*sv^2+1-tv^2);
  isb, bet := IsSquare(bq);
  tb := tv^2 - bet^2;
  isg, gam := IsSquare(tb);
  printf "  (m,n,v)=(%o,%o,%o): beta=%o  t^2-beta^2=%o square=%o", mv,nv,vv,bet,tb,isg;
  if isg then
    // recover p from beta = t(1-p^2)/(1+p^2)  =>  p^2 = (t-beta)/(t+beta)
    p2 := (tv-bet)/(tv+bet);
    isp, pv := IsSquare(p2);
    printf "  p^2=%o p_rational=%o", p2, isp;
    if isp then
      w2 := sv^2*(1+pv^2)^2 - 4*alv^2*tv^2*pv^2;
      isw, wv := IsSquare(w2);
      printf "  quartic(**)_square=%o", isw;
    end if;
  end if;
  printf "\n";
end for;
printf "CHART_VERIFY_DONE\n";
quit;
