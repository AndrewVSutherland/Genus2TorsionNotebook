// claude_ov_lane1_exact.m -- Lane 1, 2026-07-25 (resumed session).
//
// EXACT TORSION + GENERICITY for every member of the u = -1/2 contact-7
// [2,2,14] family (see code/claude_ov_lane1_family.m for the derivation).
//
// Member = a point P = k*G + eps*T on  E : y^2 = x^3 + 6x^2 - 16  (cond 288,
// E(Q) = <G=(4,-12)> x <T=(-2,0)> = Z x Z/2), via
//     m = 4/x(P),  w = y(P)*m^2/4  (so w^2 = -m^4+6m^2+4m),
//     s = (m^2+w)/d, t = (m^2-w)/d, u = -1/2,  d = (m+1)^2(m-2),
//     c4 = (G7(s)-G7(t))/(s^2-t^2), c0 = G7(s)-c4 s^2, b = c4-2, a = 9/2-c0-c4,
//     h = 1-(7/2)x+a x^2+b x^3,  f = (h^2+(x-1)^7)/x^2  (monic quintic),
// marked class D = [(1,h(1)) - infty] of order 7.
//
// For each member this script computes, WITHOUT calling TorsionSubgroup:
//   (i)   the quintic factor type  => the rational 2-rank  (deg 5: rank = #factors-1)
//   (ii)  Order(D) on the integral model (must be 7)
//         => CONTAINMENT  (Z/2)^3 x Z/7 = [2,2,14]  <=  J(Q)_tors
//   (iii) N_p = #J(F_p) = L_p(1) for every good p in [11,PMAX]
//         => J(Q)_tors divides g := gcd_p N_p   (torsion injects, p odd good)
//         g = 56  ==>  J(Q)_tors is EXACTLY [2,2,14].
//   (iv)  the RM pre-screen (squarefree core of the real-subfield disc of chi_p)
//   (v)   the strict two-prime End=Z certificate (chi irreducible, no root-power
//         degree drop for n=2..12, two such primes with linearly disjoint
//         splitting fields).
//   (vi)  a WITNESS PAIR (p,q) with gcd(N_p,N_q) = 56, recorded explicitly:
//         these two primes are what the infinite-family congruence argument
//         needs (any member congruent to this one mod p and mod q has the
//         same N_p, N_q, hence the same exact torsion).
//
// Usage: magma -b KMAX:=24 PMAX:=200 NCH:=12 code/claude_ov_lane1_exact.m

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned KMAX  then KMAX  := 24; elif Type(KMAX)  eq MonStgElt then KMAX  := StringToInteger(KMAX);  end if;
if not assigned PMAX  then PMAX  := 200; elif Type(PMAX) eq MonStgElt then PMAX  := StringToInteger(PMAX);  end if;
if not assigned NCH   then NCH   := 12; elif Type(NCH)   eq MonStgElt then NCH   := StringToInteger(NCH);   end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);
G7 := func<u0 | -(u0^5 - u0^3 - u0^2/2)/(u0+1)^2>;
E := EllipticCurve([0,6,0,0,-16]);
Gp := E![4,-12]; Tp := E![-2,0];
primes := PrimesInInterval(11, PMAX);

dir := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/lane1exact";
System("mkdir -p " cat dir);

// build the member data for (k,eps); returns ok, f, h, mm, s0, t0
function Member(k, eps)
  Pt := k*Gp + eps*Tp;
  if Pt eq E!0 then return false,_,_,_,_,_; end if;
  xc := Pt[1]/Pt[3];
  if xc eq 0 then return false,_,_,_,_,_; end if;
  mm := 4/xc; ww := (Pt[2]/Pt[3])*mm^2/4;
  d := (mm+1)^2*(mm-2);
  if d eq 0 then return false,_,_,_,_,_; end if;
  s0 := (mm^2+ww)/d; t0 := (mm^2-ww)/d; u0 := Q!(-1/2);
  vs := [s0,t0,u0];
  if #{v : v in vs} ne 3 or #{v^2 : v in vs} ne 3 then return false,_,_,_,_,_; end if;
  if (-1 in vs) or (0 in vs) then return false,_,_,_,_,_; end if;
  c4 := (G7(s0)-G7(t0))/(s0^2-t0^2); c0 := G7(s0) - c4*s0^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  num := h^2 + (x-1)^7;
  if num mod x^2 ne 0 then return false,_,_,_,_,_; end if;
  f := num div x^2;
  if Degree(f) ne 5 or Discriminant(f) eq 0 then return false,_,_,_,_,_; end if;
  return true, f, h, mm, s0, t0;
end function;

for ch in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then                                       // ---------- child ----------
    FH := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
    for k in [1..KMAX] do
      if (k mod NCH) ne ch then continue; end if;
      for eps in [0,1] do
        ok, f, h, mm, s0, t0 := Member(k, eps);
        if not ok then
          Puts(FH, Sprintf("SKIP k=%o eps=%o", k, eps)); continue;
        end if;
        // factor type of the (small) rational quintic -> 2-rank
        ftl := Sort([Degree(tt[1]) : tt in Factorization(f)]);
        nfac := #Factorization(f);
        tworank := nfac - 1;
        // integral model  x -> X/L, y -> Y/L^3
        L := LCM([Denominator(cc) : cc in Coefficients(f)]);
        F := P![ Z!(Coefficient(f,i)*L^(6-i)) : i in [0..5] ];
        C := HyperellipticCurve(F); J := Jacobian(C);
        Ppt := C ! [L*1, L^3*Evaluate(h,1)];
        D := J ! (Ppt - PointsAtInfinity(C)[1]);
        ordD := Order(D);
        // ---- N_p over good primes ----
        lc := Z!LeadingCoefficient(F);
        Nps := []; gg := 0; sig := {}; nsig := 0; certs := [];
        for p in primes do
          if lc mod p eq 0 then continue; end if;
          Fp := PolynomialRing(GF(p))!F;
          if Degree(Fp) ne 5 or not IsSquarefree(Fp) then continue; end if;
          Cp := HyperellipticCurve(Fp);
          Lp := EulerFactor(Jacobian(Cp));
          Np := Z!Evaluate(Lp, 1);
          Append(~Nps, <p, Np>);
          gg := GCD(gg, Np);
          chi := P!Reverse(Coefficients(Lp));
          dd := Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p);
          if dd ne 0 then Include(~sig, Squarefree(Z!Abs(dd))); nsig +:= 1; end if;
          if #certs lt 8 and IsIrreducible(chi) then
            K := NumberField(chi); pi := K.1; drop := false;
            for nn in [2..12] do
              if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
            end for;
            if not drop then Append(~certs, <p, chi>); end if;
          end if;
        end for;
        // witness pair with gcd exactly 56
        wp := 0; wq := 0;
        for i in [1..#Nps] do
          for j in [i+1..#Nps] do
            if GCD(Nps[i][2], Nps[j][2]) eq 56 then wp := Nps[i][1]; wq := Nps[j][1]; break; end if;
          end for;
          if wp ne 0 then break; end if;
        end for;
        // strict two-prime certificate
        cp := 0; cq := 0; cd := 0;
        for i in [1..#certs] do
          for j in [i+1..#certs] do
            d0 := Degree(SplittingField(certs[i][2]));
            d1 := Degree(SplittingField(certs[j][2]));
            dc := Degree(SplittingField(certs[i][2]*certs[j][2]));
            if dc eq d0*d1 then cp := certs[i][1]; cq := certs[j][1]; cd := dc; break; end if;
          end for;
          if cp ne 0 then break; end if;
        end for;
        ss := Sort([z : z in sig]);
        hgt := Maximum([Abs(Numerator(cc)) + Abs(Denominator(cc)) : cc in Coefficients(f)]);
        Puts(FH, Sprintf("MEMBER k=%o eps=%o m=%o s=%o t=%o log10ht=%o ftype=%o 2rank=%o ordD=%o nprimes=%o gcdNp=%o exact=%o witness=(%o,%o) ncores=%o cores=%o cert=(%o,%o,%o)",
             k, eps, mm, s0, t0, Ilog(10,Maximum(2,hgt)), ftl, tworank, ordD, #Nps, gg,
             (gg eq 56 and tworank eq 3 and ordD eq 7), wp, wq, #ss, ss[1..Min(8,#ss)], cp, cq, cd));
        Puts(FH, Sprintf("NPS k=%o eps=%o %o", k, eps, Nps));
      end for;
    end for;
    Flush(FH); delete FH;
    quit;
  end if;
end for;
WaitForAllChildren();                                     // ---------- parent ----------

lines := [];
for ch in [0..NCH-1] do
  fn := dir cat "/part" cat IntegerToString(ch) cat ".txt";
  for l in Split(Read(fn), "\n") do
    if #l gt 0 then Append(~lines, l); end if;
  end for;
end for;
nmem := 0; nexact := 0; ncert := 0; nrm := 0;
for l in lines do
  if l[1..3] eq "MEM" then
    nmem +:= 1;
    if "exact=true" in l then nexact +:= 1; end if;
    if not "cert=(0,0,0)" in l then ncert +:= 1; end if;
    printf "%o\n", l;
  end if;
end for;
printf "\n=== NPS TABLE ===\n";
for l in lines do if l[1..3] eq "NPS" then printf "%o\n", l; end if; end for;
for l in lines do if l[1..3] eq "SKI" then printf "%o\n", l; end if; end for;
printf "\nSUMMARY members=%o exact_2_2_14=%o with_strict_cert=%o KMAX=%o PMAX=%o\n", nmem, nexact, ncert, KMAX, PMAX;
printf "LANE1_EXACT_DONE\n";
quit;
