// claude_224_kernel_scan.m — [2,24] incompatibility microscope.
// For every known mech-A [2,12] point (data/claude_prod_06_224_mechA_212family*.txt),
// compute the two rational x-T coordinates of the order-12 class D on the odd quintic:
//   s1 = c*(xP - a/z)   (root rho0), s2 = c*(xP - k)  (mech-A root),
// with c = lc(f5), xP = -1/w.  D halvable requires both to be squares (plus the cubic
// norm condition).  Tally squarefree kernels: a SYSTEMATIC nonsquare pattern = the
// conjectured 2-descent incompatibility, and tells us exactly what to prove.
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

tallies1 := AssociativeArray(); tallies2 := AssociativeArray(); tallyboth := AssociativeArray();
npts := 0; ns1sq := 0; ns2sq := 0; nbothsq := 0; ncubfail := 0;
for fn in ["data/claude_prod_06_224_mechA_212family.txt", "data/claude_prod_06_224_mechA_212family_deep.txt"] do
  ok, lines := true, "";
  try lines := Read(fn); catch e ok := false; end try;
  if not ok then continue; end if;
  for ln in Split(lines, "\n") do
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
    xP := -1/w;
    rho0 := a/z;
    if Evaluate(f5, rho0) ne 0 or Evaluate(f5, k) ne 0 then continue; end if;
    npts +:= 1;
    s1 := c*(xP - rho0); s2 := c*(xP - k);
    k1 := SqfKer(s1); k2 := SqfKer(s2);
    if k1 eq 1 then ns1sq +:= 1; end if;
    if k2 eq 1 then ns2sq +:= 1; end if;
    if k1 eq 1 and k2 eq 1 then
      nbothsq +:= 1;
      // both rational coordinates square yet torsion is exactly [2,12]:
      // the cubic-part condition must be failing — record loudly
      printf "BOTHSQ t=%o/%o k=%o r=%o (cubic condition must fail)\n", tn, td, k, r;
      ncubfail +:= 1;
    end if;
    if IsDefined(tallies1, k1) then tallies1[k1] +:= 1; else tallies1[k1] := 1; end if;
    if IsDefined(tallies2, k2) then tallies2[k2] +:= 1; else tallies2[k2] := 1; end if;
    kb := SqfKer(s1*s2);
    if IsDefined(tallyboth, kb) then tallyboth[kb] +:= 1; else tallyboth[kb] := 1; end if;
  end for;
end for;
printf "points=%o  s1 square: %o  s2 square: %o  both square: %o\n", npts, ns1sq, ns2sq, nbothsq;
procedure Top(tal, name)
  arr := [ <tal[key], key> : key in Keys(tal) ];
  Sort(~arr, func<A,B | B[1]-A[1]>);
  printf "%o top kernels: %o\n", name, arr[1..Minimum(12,#arr)];
end procedure;
Top(tallies1, "s1"); Top(tallies2, "s2"); Top(tallyboth, "s1*s2");
print "ALL_DONE";
quit;
