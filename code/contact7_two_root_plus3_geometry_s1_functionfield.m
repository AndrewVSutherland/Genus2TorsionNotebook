//////////////////////////////////////////////////////////////////////
// Generic s=1 cubic-contact fiber over Q(t).  This avoids a four-variable
// elimination over Q by treating the base coordinate t as a coefficient.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=8;
elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned compute_genus then compute_genus:=false;
elif Type(compute_genus) eq MonStgElt then
    compute_genus:=compute_genus in {"true","True","1","yes"};
end if;

if not assigned p then p:=0;
elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if p eq 0 then Q:=Rationals(); else Q:=GF(p); end if;
F<t>:=FunctionField(Q);
R<U,V,M>:=PolynomialRing(F,3,"lex");
RX<x>:=PolynomialRing(R);
w:=1-t^2;
a:=F!35/8;
b:=(t^7-1+(F!7/2)*w-a*w^2)/w^3;
h:=1-(F!7/2)*x+a*x^2+b*x^3;
f:=ExactQuotient(h^2+(x-1)^7,x^2);
c:=[R!Coefficient(f,i):i in [0..5]];
assert c[1] eq 0;
A:=(M*c[6]+3*U)/2;
B:=(M*c[5]+3*(U^2+V)-A^2)/2;
E:=(M*c[4]+U^3+6*U*V-2*A*B)/2;
G2:=B^2+2*A*E-3*(U^2*V+V^2)-M*c[3];
G1:=2*B*E-3*U*V^2-M*c[2];
G0:=E^2-V^3;
I:=ideal<R|G2,G1,G0>;
// The M=0 component is the trivial identity H^2=q^3.  Saturate it.
print "S1_FUNCTION_FIELD_SAT_M_BEGIN characteristic",p;
time Is:=Saturation(I,ideal<R|M>);
time gbs:=GroebnerBasis(Is);
print "S1_FUNCTION_FIELD_SAT_M_DONE",#gbs,
      [<TotalDegree(g),Degree(g,U),Degree(g,V),Degree(g,M),#Terms(g)>:g in gbs];
univ:=[g:g in gbs|Degree(g,U) eq 0 and Degree(g,V) eq 0];
print "M_UNIVARIATES",#univ;
for g in univ do
    print " M_POLY",<Degree(g,M),#Terms(g)>;
    print " M_FACTORS",[<Degree(fe[1],M),fe[2]>:fe in Factorization(g)];
    PM<Y>:=PolynomialRing(F);
    pm:=&+[PM!(Coefficient(g,M,i))*Y^i:i in [0..Degree(g,M)]];
    if #Factorization(pm) eq 1 and Factorization(pm)[1][2] eq 1 then
        print " QUOTIENT_FUNCTION_FIELD_BEGIN";
        n:=Degree(pm);
        D:=LCM([Denominator(Coefficient(pm,i)):i in [0..n]]);
        pint:=&+[PM!(Coefficient(pm,i)*D^(n-i))*Y^i:i in [0..n]];
        assert IsMonic(pint) and &and[Denominator(cc) eq 1:cc in Coefficients(pint)];
        print " INTEGRAL_SCALE_DEGREE",Degree(D);
        K40<zroot>:=FunctionField(pint);
        mroot:=zroot/K40!D;
        PK40<Lvar>:=PolynomialRing(K40);
        sq:=Lvar^2-mroot;
        print " EXACT_SQUARE_COVER_IRREDUCIBLE",IsIrreducible(sq);
        if compute_genus then
            time g40,ci40:=Genus(K40 : Al:="Montes", IsExact:=true);
            print " QUOTIENT_FUNCTION_FIELD_GENUS",g40,"constant_index",ci40;
        end if;
        if compute_genus and IsIrreducible(sq) then
            K80<lroot>:=FunctionField(sq);
            time g80,ci80:=Genus(K80 : Al:="Montes", IsExact:=true);
            print " EXACT_FUNCTION_FIELD_GENUS",g80,"constant_index",ci80;
        end if;
    end if;
end for;
quit;
