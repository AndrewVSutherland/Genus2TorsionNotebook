//////////////////////////////////////////////////////////////////////
// Rational parametrization of the square cover of the exact
// shared-Weierstrass point-order-10 locus on compact M(12).
//
// In the notation of m12_shared_weierstrass_order10_geometry.m, the
// constant coefficient equation is
//
//   s*R^5=-4*b^3*(b-1)^2.
//
// Hence s is a square iff -b*R is a square.  Put
//
//   R=-z^2/b, a=2*b^4*(b-1)/z^5, s=a^2.
//
// The two highest remaining equations solve linearly for C=a*c and
// D=a*d.  The last two equations below cut out the actual rational
// point-order-10 curve in (b,w,z).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","resultant"} then
    error "mode must be summary or resultant";
end if;
if not assigned MemGB then MemGB := 3;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); Zint := Integers();
R<b,w,z> := PolynomialRing(Q,3,"grevlex");
K := FieldOfFractions(R);

g0 := -4*b^5+8*b^4-4*b^3;
g1 := 4*b^6*w^2-20*b^5*w^2+24*b^5*w+41*b^4*w^2
      -72*b^4*w-44*b^3*w^2+4*b^4+78*b^3*w+26*b^2*w^2
      -12*b^3-36*b^2*w-8*b*w^2+9*b^2+6*b*w+w^2;
g2 := -8*b^5*w^2+24*b^4*w^2+16*b^4*w-26*b^3*w^2
      -32*b^3*w+12*b^2*w^2+12*b^2*w-2*b*w^2+4*b*w
      -6*b-2*w;
g3 := 4*b^4*w^2-4*b^3*w^2+24*b^3*w+b^2*w^2-24*b^2*w
      +4*b^2+6*b*w+4*b+1;
g4 := -4*b;

Rp := -z^2/b;
a := 2*b^4*(b-1)/z^5;
s := a^2;
C := (g4-5*s*Rp)/2;
D := (s*g3-C^2+10*s^2*Rp^2)/(2*s);
E2rat := s*g2-2*C*D-10*s^2*Rp^3;
E1rat := s*g1-D^2+5*s^2*Rp^4;

function PrimitiveNumerator(f)
    n := R!Numerator(f);
    den := LCM([Denominator(c) : c in Coefficients(n)]);
    nz := [Zint!(den*c) : c in Coefficients(n)];
    cont := GCD(nz); if cont eq 0 then cont:=1; end if;
    return R!((Q!den/Q!Abs(cont))*n);
end function;

function StripBoundary(f)
    changed := true;
    while changed do
        changed := false;
        for h in [b,b-1,z] do
            if IsDivisibleBy(f,h) then
                f := ExactQuotient(f,h); changed := true;
            end if;
        end for;
    end while;
    return f;
end function;

N2raw := PrimitiveNumerator(E2rat);
N1raw := PrimitiveNumerator(E1rat);
N2 := StripBoundary(N2raw);
N1 := StripBoundary(N1raw);

PX<X> := PolynomialRing(K);
h := a*X^2+(C/a)*X+(D/a);
G := g0+g1*X+g2*X^2+g3*X^3+g4*X^4;
rec := G-X*h^2+s*(X-Rp)^5;
assert Coefficient(rec,5) eq 0;
assert &and[Coefficient(rec,i) eq 0 : i in [0,3,4]];
assert IsDivisibleBy(Numerator(Coefficient(rec,2)),N2);
assert IsDivisibleBy(Numerator(Coefficient(rec,1)),N1);

print "M12_SHARED_W_PARAM_SELF_TEST_PASS";
print "PARAMETERS","R=-z^2/b","a=2*b^4*(b-1)/z^5";
print "N2_SHAPE",<TotalDegree(N2),Degree(N2,b),Degree(N2,w),Degree(N2,z),#Terms(N2)>;
print "N1_SHAPE",<TotalDegree(N1),Degree(N1,b),Degree(N1,w),Degree(N1,z),#Terms(N1)>;
print "N2_FACTORIZATION";
for fe in Factorization(N2) do
    print fe[2],<TotalDegree(fe[1]),#Terms(fe[1])>,fe[1];
end for;
print "N1_FACTORIZATION";
for fe in Factorization(N1) do
    print fe[2],<TotalDegree(fe[1]),#Terms(fe[1])>,fe[1];
end for;
if mode eq "summary" then quit; end if;

print "RESULTANT_W_BEGIN";
time Res := Resultant(N2,N1,w);
Res := PrimitiveNumerator(Res);
print "RESULTANT_W_SHAPE",<TotalDegree(Res),Degree(Res,b),Degree(Res,z),#Terms(Res)>;
print "RESULTANT_W_FACTORIZATION";
for fe in Factorization(Res) do
    print fe[2],<TotalDegree(fe[1]),Degree(fe[1],b),Degree(fe[1],z),#Terms(fe[1])>;
    if #Terms(fe[1]) le 120 then print fe[1]; end if;
    if Degree(fe[1],z) gt 1 then
        lcz := Coefficient(fe[1],z,Degree(fe[1],z));
        print "NONBOUNDARY_Z_LEADING_FACTORS",Factorization(lcz);
        print "NONBOUNDARY_Z_LEADING_SCALAR",lcz/(b-1/2)^6;
        print "NONBOUNDARY_Z_CONSTANT_FACTORS",Factorization(Coefficient(fe[1],z,0));
        print "NONBOUNDARY_Z_CONSTANT_SCALAR",Coefficient(fe[1],z,0)/(b^42*(b-1)^22);
        assert &and[IsEven(Degree(mon,z)) : mon in Monomials(fe[1])];
        S<bb,t> := PolynomialRing(Q,2);
        phi := hom<R -> S | bb,S!0,S!0>;
        evenModel := &+[phi(Coefficient(fe[1],z,2*i))*t^i : i in [0..Degree(fe[1],z) div 2]];
        print "EVEN_QUOTIENT_SHAPE",<TotalDegree(evenModel),Degree(evenModel,bb),Degree(evenModel,t),#Terms(evenModel)>;
        print "EVEN_QUOTIENT_FACTORS",[<TotalDegree(ff[1]),Degree(ff[1],bb),Degree(ff[1],t),#Terms(ff[1]),ff[2]> : ff in Factorization(evenModel)];
    end if;
end for;
quit;
