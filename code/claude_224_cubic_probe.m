// claude_224_cubic_probe.m — identify the mechanism of the [2,24] cubic-coordinate
// obstruction: at mech-A points where both rational x-T coordinates are squares,
// probe epsilon = c*(xP - theta) in the cubic field L = Q[theta]:
//   - IsSquare(epsilon) (expected false),
//   - Norm(epsilon) square-class (expected square: product relation),
//   - REAL signs of epsilon at the real embeddings of L (a constant sign pattern
//     = provable real obstruction),
//   - local nonsquare places (odd p dividing the relevant discriminants, and 2).
// Also: verify symbolically that s1 = c*(xP - rho0) is a square identically.
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals(); P<X> := PolynomialRing(Q);

// ---- symbolic s1-identity over Q(z, r) ----
Pzr<zz, rr> := PolynomialRing(Q, 2);
Fzr := FieldOfFractions(Pzr);
PF<XX> := PolynomialRing(Fzr);
az := (1-zz^2)/(4*(rr+1)); wz := 2*(rr+1)/(1+zz);
Tz := az*XX^2 - XX + rr; hz := (XX-rr)*(Tz+1);
Wz := hz^2 + 4*az*XX^2*Tz*(Tz+1);
f5z := PF!0;
for i in [0..Degree(Wz)] do for j in [0..i] do
  f5z +:= Coefficient(Wz, i)*Binomial(i,j)*wz^(i-j)*XX^(6-j);
end for; end for;
cz := LeadingCoefficient(f5z);
s1z := cz*(-1/wz - az/zz);
nu := Numerator(s1z)*Denominator(s1z);
fac := Factorization(Numerator(nu) * Denominator(nu));
pr := &*[ Parent(nu) | t[1]^t[2] : t in fac ];
unit := nu div pr;
odd := &*[ Parent(nu) | t[1]^(t[2] mod 2) : t in fac ];
printf "s1 symbolic: oddpart=%o unit=%o => identically square: %o\n",
  odd, unit, (odd eq 1) and IsSquare(Q!unit);

// ---- probe the BOTHSQ points ----
function SqfKer(v)
  if v eq 0 then return 0; end if;
  n := Numerator(v)*Denominator(v);
  s := Sign(n); n := AbsoluteValue(n);
  k := 1;
  for pf in Factorization(n) do if IsOdd(pf[2]) then k *:= pf[1]; end if; end for;
  return s*k;
end function;

nprobe := 0;
for fn in ["data/claude_prod_06_224_mechA_212family.txt", "data/claude_prod_06_224_mechA_212family_deep.txt"] do
  ok, lines := true, "";
  try lines := Read(fn); catch e ok := false; end try;
  if not ok then continue; end if;
  for ln in Split(lines, "\n") do
    if nprobe ge 40 then break fn; end if;
    if #ln eq 0 or ln[1] eq "#" then continue; end if;
    w6 := Split(ln, " ");
    if #w6 lt 6 then continue; end if;
    tn := StringToInteger(w6[1]); td := StringToInteger(w6[2]);
    kn := StringToInteger(w6[3]); kd := StringToInteger(w6[4]);
    rn := StringToInteger(w6[5]); rd := StringToInteger(w6[6]);
    z := Q!(tn^2+td^2)/(2*tn*td); r := Q!rn/rd; k := Q!kn/kd;
    if r eq -1 or r eq 0 or z^2 eq 1 then continue; end if;
    a := (1-z^2)/(4*(r+1)); w := 2*(r+1)/(1+z);
    T := a*X^2 - X + r; h := (X-r)*(T+1);
    W := h^2 + 4*a*X^2*T*(T+1);
    f5 := P!0;
    for i in [0..Degree(W)] do for j in [0..i] do
      f5 +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
    end for; end for;
    if Degree(f5) ne 5 or Discriminant(f5) eq 0 then continue; end if;
    c := LeadingCoefficient(f5);
    xP := -1/w; rho0 := a/z;
    s2 := c*(xP - k);
    if SqfKer(s2) ne 1 then continue; end if;
    // cubic factor
    g3 := f5 div ((X - rho0)*(X - k)*c);
    g3 := g3 / LeadingCoefficient(g3);
    if Degree(g3) ne 3 or not IsIrreducible(g3) then
      printf "SPLITCUBIC t=%o/%o k=%o r=%o g3 %o\n", tn, td, k, r, Factorization(g3);
      continue;
    end if;
    nprobe +:= 1;
    L<th> := NumberField(g3);
    eps := c*(xP - th);
    issq := IsSquare(eps);
    nrm := SqfKer(Q!Norm(eps));
    sig := Signature(L);
    realsigns := [ Sign(Evaluate(c*(xP - cj), 0)) : cj in [] ];  // placeholder
    cjs := Conjugates(eps);
    rsig := [ Sign(Re(cjs[i])) : i in [1..sig] ];
    printf "PROBE t=%o/%o k=%o r=%o : eps square %o, Norm kernel %o, sig(L)=(%o), real signs %o\n",
      tn, td, k, r, issq, nrm, sig, rsig;
  end for;
end for;
printf "probed %o both-square points\n", nprobe;
print "ALL_DONE";
quit;
