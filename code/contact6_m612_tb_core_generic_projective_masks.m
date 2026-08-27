//////////////////////////////////////////////////////////////////////
// Conservative projective masks for the fixed (s,v) slices of the
// T_B-halving + [3,3] fiber product.
//
// The crucial order of operations is:
//   (1) work over Q(s,v),
//   (2) impose K=m^2 and M=L^2,
//   (3) saturate the generic fiber by the open-chart boundary,
//   (4) take its projective closure in (r,m,L,U),
//   (5) specialize the resulting equations modulo p.
//
// Thus boundary limits of the saturated characteristic-zero family are
// retained.  A residue at which any coefficient denominator vanishes is
// marked exceptional and retained automatically.  When these masks are
// applied to rational (s,v), a denominator divisible by p must likewise be
// retained as base-projective infinity.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned primes then
    prime_list:=[5,7,11];
elif Type(primes) eq MonStgElt then
    prime_list:=[StringToInteger(t):t in Split(primes,",")|#t gt 0];
else
    prime_list:=primes;
end if;
if not assigned MemGB then MemGB:=8;
elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q:=Rationals();
B<s,v>:=PolynomialRing(Q,2);
K:=FieldOfFractions(B);
A<r,m,L,U>:=PolynomialRing(K,4,"grevlex");
FA:=FieldOfFractions(A);
PA<a>:=PolynomialRing(FA);
sF:=FA!(K!s); vF:=FA!(K!v);
rF:=FA!r; mF:=FA!m; LF:=FA!L; UF:=FA!U;
KF:=mF^2; MF:=LF^2; bF:=2*sF^2-3;

A3:=a-3+4*sF^2*rF-2*KF;
H1:=sF*((a-3)*rF^2+4*rF-KF*(a+3))-rF*A3;
H2:=8*sF^2*(2*sF^2*rF^2+2*(a-3)*rF+2-KF*(2*sF^2-6))
    -A3^2-32*sF^3*rF;
D:=Coefficient(H1,1);
N:=-Coefficient(H1,0);

function SubNum(poly)
    d:=Degree(poly);
    val:=FA!0;
    for i in [0..d] do
        val+:=Coefficient(poly,i)*N^i*D^(d-i);
    end for;
    return A!Numerator(val);
end function;

G0:=SubNum(H2);
while IsDivisibleBy(G0,m) do G0:=ExactQuotient(G0,m); end while;

c1:=2*a+6;
c2:=a^2+2*bF-15;
c3:=2*a*bF+22;
c4:=2*a+bF^2-15;
c5:=2*bF+6;
B3:=c5*MF+3*UF;
Delta3:=4*c4*MF+12*(UF^2+vF^2)-B3^2;
F3:=B3*Delta3+16*vF^3-8*c3*MF-8*UF^3-48*UF*vF^2;
F2:=Delta3^2+64*B3*vF^3-64*c2*MF-192*(UF^2*vF^2+vF^4);
F1:=Delta3*vF^3-4*c1*MF-12*UF*vF^4;
G1:=SubNum(F1); G2:=SubNum(F2); G3:=SubNum(F3);

DA:=A!Numerator(D); NA:=A!Numerator(N);
boundary:=r*m*L*DA*(U^2-4*(K!v)^2)*(NA+(2*(K!s)^2-1)*DA);
I:=ideal<A|G0,G1,G2,G3>;

print "CONTACT6_M612_TB_CORE_GENERIC_PROJECTIVE_MASKS";
print "generic_shapes",[<TotalDegree(g),#Terms(g)>:g in [G0,G1,G2,G3]];
print "generic_saturation_begin";
time Isat:=Saturation(I,ideal<A|boundary>);
print "generic_saturation_done","dimension",Dimension(Isat),
      "basis_size",#Basis(Isat);
X:=Scheme(AffineSpace(A),Basis(Isat));
print "generic_projective_closure_begin";
time Xbar:=ProjectiveClosure(X);
Jbar:=Ideal(Xbar);
eqns:=Basis(Jbar);
print "generic_projective_closure_done","dimension",Dimension(Xbar),
      "equation_count",#eqns,
      "equation_shapes",[<TotalDegree(g),#Terms(g)>:g in eqns];

// Evaluate a Q[s,v] polynomial modulo p without ever coercing a rational
// coefficient whose denominator is divisible by p.
function EvalBaseModP(poly,s0,v0,F,p)
    ans:=F!0;
    mons:=Monomials(poly); coeffs:=Coefficients(poly);
    for i in [1..#mons] do
        q:=Q!coeffs[i];
        if Denominator(q) mod p eq 0 then return false,F!0; end if;
        cq:=(F!Numerator(q))/(F!Denominator(q));
        ee:=Exponents(mons[i]);
        ans+:=cq*s0^ee[1]*v0^ee[2];
    end for;
    return true,ans;
end function;

function SpecializeEquation(g,s0,v0,F,p,Hp)
    ans:=Hp!0;
    mons:=Monomials(g); coeffs:=Coefficients(g);
    for i in [1..#mons] do
        c:=K!coeffs[i];
        okn,cn:=EvalBaseModP(Numerator(c),s0,v0,F,p);
        okd,cd:=EvalBaseModP(Denominator(c),s0,v0,F,p);
        if not okn or not okd or cd eq 0 then return false,Hp!0; end if;
        ee:=Exponents(mons[i]);
        mon:=Hp!1;
        for j in [1..#ee] do mon*:=Hp.j^ee[j]; end for;
        ans+:=(cn/cd)*mon;
    end for;
    return true,ans;
end function;

for p in prime_list do
    F:=GF(p); Hp:=PolynomialRing(F,5,"grevlex");
    allowed:=[]; exceptional:=[]; excluded:=[];
    t0:=Cputime();
    for s0 in F do
        for v0 in F do
            specs:=[]; ok:=true;
            for g in eqns do
                good,gg:=SpecializeEquation(g,s0,v0,F,p,Hp);
                if not good then ok:=false; break; end if;
                if gg ne 0 then Append(~specs,gg); end if;
            end for;
            key:=<Integers()!s0,Integers()!v0>;
            if not ok then
                Append(~exceptional,key); Append(~allowed,key);
                continue;
            end if;
            P4:=ProjectiveSpace(F,4);
            CR:=CoordinateRing(P4);
            Y:=Scheme(P4,[Evaluate(g,[CR.i:i in [1..5]]):g in specs]);
            if #Points(Y) gt 0 then Append(~allowed,key);
            else Append(~excluded,key); end if;
        end for;
    end for;
    print "p",p,"allowed",#allowed,"exceptional",#exceptional,
          "excluded",#excluded,"seconds",Cputime(t0);
    print "ALLOWED",allowed;
    print "EXCEPTIONAL",exceptional;
    print "EXCLUDED",excluded;
end for;
quit;
