//////////////////////////////////////////////////////////////////////
// Reciprocal-t model and bounded rational-slice search for the exact
// point-order-10 (= repeated-q/shared-Weierstrass) locus on compact M(12).
//
// This performs no large Groebner computation.  The w-resultant takes
// less than a second and isolates one irreducible candidate R(b,t) of
// bidegree (22,58).  In search mode, every rational t of height <= H is
// specialized and all rational b-roots are computed exactly; b is not
// height-bounded.
//
// Examples:
//   magma -b mode:="summary" code/m12_repeated_q_reciprocal_search.m
//   magma -b mode:="search" H:=50 code/m12_repeated_q_reciprocal_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
require mode in {"summary","search"}: "mode must be summary or search";
if not assigned H then H := 50;
elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MemGB then MemGB := 3;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
R<b,w,t> := PolynomialRing(Q,3,"grevlex");
K := FieldOfFractions(R);
P<z> := PolynomialRing(K);

d := 2*b-1;
Hn := z-b+w*(d+b*(z-b));
G := P!(d^2*z*Hn^2+4*b*(z+b-1)^2*(w*z*d^2-(z-b)^2));
g := [K!Coefficient(G,i):i in [0..4]];

r := -b/t^2;
c2 := 2*(b-1)*t^5/b;
kappa := -c2^2;
c1 := (g[5]+5*kappa*r)/(2*c2);
c0 := (g[3]+10*kappa*r^3)/(2*c1);

N1 := R!Numerator(g[2]-c0^2-5*kappa*r^4);
N3 := R!Numerator(g[4]-c1^2-2*c0*c2-10*kappa*r^2);

print "M12_REPEAT_RECIPROCAL_FORMULAS_PASS";
print "N1_SHAPE",<TotalDegree(N1),Degree(N1,b),Degree(N1,w),
                  Degree(N1,t),#Terms(N1)>;
print "N3_SHAPE",<TotalDegree(N3),Degree(N3,b),Degree(N3,w),
                  Degree(N3,t),#Terms(N3)>;

time resultant := Resultant(N1,N3,w);
fac := Factorization(resultant);
print "RESULTANT_FACTORS",
      [<TotalDegree(a[1]),Degree(a[1],b),Degree(a[1],t),
        #Terms(a[1]),a[2]>:a in fac];

quadratic := [a[1]:a in fac|Degree(a[1],b) eq 2][1];
candidate := [a[1]:a in fac|Degree(a[1],b) eq 22][1];
expectedQuadratic := (5*t^8-1)*b^2-10*t^8*b+5*t^8;
assert expectedQuadratic eq 5*quadratic;
qdisc := Coefficient(quadratic,b,1)^2-4*Coefficient(quadratic,b,2)*Coefficient(quadratic,b,0); assert qdisc eq 4*t^8/5;
assert &and[IsEven(Degree(m,t)):m in Monomials(candidate)];
assert #Factorization(candidate) eq 1;
print "POINTLESS_QUADRATIC",
      "(5*t^8-1)*b^2-10*t^8*b+5*t^8; discriminant=20*t^8";
print "CANDIDATE_SHAPE",<TotalDegree(candidate),Degree(candidate,b),
                         Degree(candidate,t),#Terms(candidate)>;

if mode eq "summary" then quit; end if;

QB<X> := PolynomialRing(Q);
values := [];
for den in [1..H] do
    for num in [-H..H] do
        if num ne 0 and GCD(Abs(num),den) eq 1 then
            Append(~values,Q!num/den);
        end if;
    end for;
end for;

function SpecializeInB(f,tv)
    ans := QB!0;
    for i in [0..Degree(f,b)] do
        ans +:= Evaluate(Coefficient(f,b,i),[Q!0,Q!0,tv])*X^i;
    end for;
    return ans;
end function;

function SpecializeInW(f,bv,tv)
    ans := QB!0;
    for i in [0..Degree(f,w)] do
        ans +:= Evaluate(Coefficient(f,w,i),[bv,Q!0,tv])*X^i;
    end for;
    return ans;
end function;

hits := [];
checked := 0;
timer := Cputime();
for tv in values do
    fb := SpecializeInB(candidate,tv);
    if fb eq 0 then
        print "IDENTICALLY_ZERO_SLICE",tv;
        continue;
    end if;
    for rb in Roots(fb) do
        bv := rb[1];
        if bv in {Q!0,Q!1,Q!1/2} then continue; end if;
        fw := GCD(SpecializeInW(N1,bv,tv),SpecializeInW(N3,bv,tv));
        wr := Roots(fw);
        Append(~hits,<tv,bv,rb[2],wr>);
        print "RATIONAL_RESULTANT_HIT","t",tv,"b",bv,"w_roots",wr;
    end for;
    checked +:= 1;
    if checked mod 5000 eq 0 then
        print "PROGRESS",checked,#values,"cpu",Cputime(timer),"hits",#hits;
    end if;
end for;

print "RECIPROCAL_SLICE_SEARCH_DONE","H",H,"slices",checked,
      "cpu",Cputime(timer),"hits",hits;
quit;
