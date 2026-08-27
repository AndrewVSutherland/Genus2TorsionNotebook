//////////////////////////////////////////////////////////////////////
// Direct exact a=0 fiber, which is e=infinity and therefore cannot be
// evaluated through the rational recovery formula in e.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Q:=Rationals(); Qz<z>:=PolynomialRing(Q);
R<M,U,v>:=PolynomialRing(Q,3,"grevlex");
a:=Q!0;
N:=(3*U+6*M)/2;
Rc:=(3*U^2+3*v^2+(2*a-15)*M-N^2)/2;
Fs:=[2*v^3+2*N*Rc-U^3-6*U*v^2-22*M,
     Rc^2+2*N*v^3-3*U^2*v^2-3*v^4-(a^2-15)*M,
     2*Rc*v^3-3*U*v^4-(2*a+6)*M];
I:=Saturation(ideal<R|Fs>,ideal<R|M>);
print "A0_FIBER_DIMENSION",Dimension(I),
      "ALGEBRA_DIMENSION",Dimension(quo<R|I>);
GL:=GroebnerBasis(ChangeOrder(I,"lex"));
univ:=[g:g in GL|Degree(g,1) eq 0 and Degree(g,2) eq 0 and
                   Degree(g,3) gt 0];
assert #univ eq 1;
rv:=Qz!UnivariatePolynomial(univ[1]); rv/:=LeadingCoefficient(rv);
print "A0_V_FACTORS",[<Degree(h[1]),h[2],h[1]>:h in Factorization(rv)];

// Coefficientwise e=infinity (a=0) limit of the exact degree-12 factor.
f12inf:=z^12+6*z^10+10*z^9-15*z^8-18*z^7-18*z^5-15*z^4
        +10*z^3+6*z^2+1;
fac12inf:=Factorization(f12inf);
print "A0_LIMITING_DEGREE12_FACTORS",
      [<Degree(h[1]),h[2],h[1]>:h in fac12inf];
assert [<Degree(h[1]),h[2]>:h in fac12inf] eq [<1,4>,<8,1>];

gm:=[g:g in GL|Degree(g,1) eq 1 and Degree(g,2) eq 0];
gu:=[g:g in GL|Degree(g,1) eq 0 and Degree(g,2) eq 1];
print "A0_M_LINEAR_RELATIONS",#gm,"U_LINEAR_RELATIONS",#gu;
if #gm gt 0 and #gu gt 0 then
    cm:=Q!Coefficient(gm[1],1,1); cu:=Q!Coefficient(gu[1],2,1);
    m0:=(-Qz!UnivariatePolynomial(Evaluate(gm[1],[R!0,R!0,v]))/cm)
         mod rv;
    u0:=(-Qz!UnivariatePolynomial(Evaluate(gu[1],[R!0,R!0,v]))/cu)
         mod rv;
    hits:=[];
    for rt in Roots(rv) do
        vv:=rt[1]; mv:=Evaluate(m0,vv); uv:=Evaluate(u0,vv);
        sq,Lv:=IsSquare(mv);
        if sq then Append(~hits,<vv,mv,Lv,uv,
            mv ne 0 and vv ne 0 and uv^2-4*vv^2 ne 0>); end if;
    end for;
    print "A0_RATIONAL_SIGNED_HITS",hits;
end if;

// The degree-12 support factor specializes to (v+1)^4 times the displayed
// degree-8 factor.  Resolve its only rational support value v=-1 directly.
Jminus:=I+ideal<R|v+1>;
Gminus:=GroebnerBasis(ChangeOrder(Jminus,"lex"));
print "A0_V_MINUS1_LEX_BASIS",Gminus;
uonly:=[g:g in Gminus|Degree(g,1) eq 0 and Degree(g,2) gt 0 and
                     Degree(g,3) eq 0];
assert #uonly eq 1;
upol:=Qz!UnivariatePolynomial(uonly[1]);
print "A0_V_MINUS1_U_FACTORS",Factorization(upol),
      "RATIONAL_U_ROOTS",Roots(upol);
assert #Roots(upol) eq 0;
quit;
