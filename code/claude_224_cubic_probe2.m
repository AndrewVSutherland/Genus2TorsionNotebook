// claude_224_cubic_probe2.m — pinpoint the finite-place mechanism of the [2,24]
// cubic obstruction: for both-square mech-A points, factor the ideal (eps) in O_L,
// record primes with odd valuation (residue characteristics), and test 2-adic
// squareness at the primes above 2.  A uniform pattern identifies the mechanism.
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals(); P<X> := PolynomialRing(Q);
function SqfKer(v)
  if v eq 0 then return 0; end if;
  n := Numerator(v)*Denominator(v);
  s := Sign(n); n := AbsoluteValue(n);
  k := 1;
  for pf in Factorization(n) do if IsOdd(pf[2]) then k *:= pf[1]; end if; end for;
  return s*k;
end function;

nprobe := 0;
for fn in ["data/claude_prod_06_224_mechA_212family.txt"] do
  lines := Read(fn);
  for ln in Split(lines, "\n") do
    if nprobe ge 12 then break fn; end if;
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
    if SqfKer(c*(xP - k)) ne 1 then continue; end if;
    g3 := f5 div ((X - rho0)*(X - k)*c);
    g3 := g3 / LeadingCoefficient(g3);
    if Degree(g3) ne 3 or not IsIrreducible(g3) then continue; end if;
    nprobe +:= 1;
    L<th> := NumberField(g3);
    OL := MaximalOrder(L);
    eps := c*(xP - th);
    fac := Factorization(eps*OL);
    oddpr := [ <Norm(t[1]), MinimalInteger(t[1]), t[2]> : t in fac | IsOdd(t[2]) ];
    // 2-adic squareness at primes above 2
    twos := [];
    for pr in Decomposition(OL, 2) do
      P2 := pr[1];
      loc := true;
      try
        C2, mp2 := Completion(L, P2 : Precision := 60);
        loc := IsSquare(mp2(eps));
      catch e loc := true; end try;
      Append(~twos, loc);
    end for;
    printf "PROBE2 t=%o/%o k=%o : odd-val primes (norm,p,val)=%o ; square at 2-adic places %o\n",
      tn, td, k, oddpr, twos;
  end for;
end for;
printf "probed %o\n", nprobe;
print "ALL_DONE";
quit;
