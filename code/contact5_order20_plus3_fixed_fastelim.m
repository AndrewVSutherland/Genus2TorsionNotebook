//////////////////////////////////////////////////////////////////////
// Fast elimination in the fixed-Weierstrass square-cover chart.
// This chart already imposes M=L^2.  The two equations F4,F2 are
// quadratic in a; recover a=N/D exactly and substitute into F3.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=8; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q:=Rationals();
C<s,r,j>:=PolynomialRing(Q,3,"grevlex");
PA<a>:=PolynomialRing(C);

F4 := r^5*j^4 - 4*s^2*r^3*j^2 - 6*r^3*a*j^2 + 20*r^3*j^2
      -4*s^2*j^2 -12*r^3+12*r^2*a-3*r*a^2
      +16*s*j^2-16*j^2;
F3 := 2*s^3*r^3*j^2+s^2*r^2*j^4-8*s^2*r^3*j^2
      -3*r^4*a*j^2-4*s*r^2*j^4-2*s*r^3*j^2+4*r^2*j^4
      +20*r^3*j^2-3*s^2*a*j^2+4*r^4-12*r^3*a+9*r^2*a^2
      -2*r*a^3+12*s*a*j^2-12*a*j^2;
F2 := -s^4*r^4*j^2+8*s^3*r^4*j^2-4*r^7*j^2-14*s^2*r^4*j^2
      +s^4*j^4-16*s*r^4*j^2-6*s^2*r^2*a*j^2-8*s^3*j^4
      -12*r^6+12*r^5*a-3*r^4*a^2+39*r^4*j^2
      +24*s*r^2*a*j^2+24*s^2*j^4-24*r^2*a*j^2
      -32*s*j^4+16*j^4;
assert <Degree(F4),Degree(F2),Degree(F3)> eq <2,2,3>;

a2:=Coefficient(F4,2); a1:=Coefficient(F4,1); a0:=Coefficient(F4,0);
b2:=Coefficient(F2,2); b1:=Coefficient(F2,1); b0:=Coefficient(F2,0);
D:=b2*a1-a2*b1;
N:=a2*b0-b2*a0;
X:=a1*b0-a0*b1;
Rcore:=N^2+D*X;
function ClearEval(f,n,d)
    return &+[Coefficient(f,i)*n^i*d^(Degree(f)-i):i in [0..Degree(f)]];
end function;
Q4:=ClearEval(F4,N,D);
Q2:=ClearEval(F2,N,D);
Q3:=ClearEval(F3,N,D);
assert Q4 eq a2*Rcore;
assert Q2 eq b2*Rcore;

print "FIXED_FAST_ELIMINATION";
print "RECOVERY","D",<TotalDegree(D),#Terms(D)>,
      "N",<TotalDegree(N),#Terms(N)>,
      "R",<TotalDegree(Rcore),#Terms(Rcore),Degree(Rcore,r),Degree(Rcore,j)>,
      "Q3",<TotalDegree(Q3),#Terms(Q3),Degree(Q3,r),Degree(Q3,j)>;

common:=GCD(Rcore,Q3);
print "COMMON_FACTOR",<TotalDegree(common),#Terms(common)>,
      "factorization",Factorization(common);
Ropen:=ExactQuotient(Rcore,common);
Q3open:=ExactQuotient(Q3,common);
print "OPEN_RECOVERY_EQUATIONS",
      <TotalDegree(Ropen),#Terms(Ropen),Degree(Ropen,r),Degree(Ropen,j)>,
      <TotalDegree(Q3open),#Terms(Q3open),Degree(Q3open,r),Degree(Q3open,j)>;

print "RESULTANT_R_BEGIN";
time ResR:=Resultant(Ropen,Q3open,r);
print "RESULTANT_R_DONE",<TotalDegree(ResR),#Terms(ResR)>;
print "RESULTANT_J_BEGIN";
time ResJ:=Resultant(Ropen,Q3open,j);
print "RESULTANT_J_DONE",<TotalDegree(ResJ),#Terms(ResJ)>;

procedure Strip(~f,h,label)
    n:=0; while IsDivisibleBy(f,h) do f:=ExactQuotient(f,h); n+:=1; end while;
    if n gt 0 then print "STRIPPED",label,n; end if;
end procedure;

// Both plane projections retain the family-boundary factors for explicit
// classification; strip only the obvious zero square parameters.
CR<sR,jR>:=PolynomialRing(Q,2,"grevlex");
toR:=hom<C->CR|sR,CR!0,jR>;
PR:=CR!toR(ResR);
Strip(~PR,jR,"j=0 in (s,j)");
Strip(~PR,sR-1,"s=1 in (s,j)");
Strip(~PR,sR-2,"s=2 in (s,j)");
Strip(~PR,8*sR^3-59*sR^2-18*sR+197,"Ds in (s,j)");
print "OPEN_SJ",<TotalDegree(PR),#Terms(PR)>;
time facR:=Factorization(PR);
print "SJ_FACTORS",#facR;
for fe in facR do
 print " SJ_FACTOR",fe[2],TotalDegree(fe[1]),Degree(fe[1],sR),
       Degree(fe[1],jR),#Terms(fe[1]);
 if TotalDegree(fe[1]) le 5 then print "  SJ_SMALL",fe[1]; end if;
end for;

CJ<sJ,rJ>:=PolynomialRing(Q,2,"grevlex");
toJ:=hom<C->CJ|sJ,rJ,CJ!0>;
PJ:=CJ!toJ(ResJ);
Strip(~PJ,rJ,"r=0 in (s,r)");
Strip(~PJ,sJ-1,"s=1 in (s,r)");
Strip(~PJ,sJ-2,"s=2 in (s,r)");
Strip(~PJ,8*sJ^3-59*sJ^2-18*sJ+197,"Ds in (s,r)");
print "OPEN_SR",<TotalDegree(PJ),#Terms(PJ)>;
time facJ:=Factorization(PJ);
print "SR_FACTORS",#facJ;
mainSR:=CJ!0;
for fe in facJ do
 print " SR_FACTOR",fe[2],TotalDegree(fe[1]),Degree(fe[1],sJ),
       Degree(fe[1],rJ),#Terms(fe[1]);
 if TotalDegree(fe[1]) le 5 then print "  SR_SMALL",fe[1]; end if;
 if TotalDegree(fe[1]) eq 44 then mainSR:=fe[1]; end if;
end for;

//////////////////////////////////////////////////////////////////////
// Recover k=j^2 on the degree-44 quotient.
//////////////////////////////////////////////////////////////////////

assert mainSR ne 0;
PK<kvar>:=PolynomialRing(CJ);
toJ0:=hom<C->CJ|sJ,rJ,CJ!0>;
function EvenJToK(f)
    ans:=PK!0;
    dj:=Degree(f,3);
    for e in [0..dj] do
        ce:=Coefficient(f,3,e);
        if IsOdd(e) then assert ce eq 0; continue; end if;
        ans +:= toJ0(ce)*kvar^(e div 2);
    end for;
    return ans;
end function;
RK:=EvenJToK(Ropen);
QK:=EvenJToK(Q3open);
print "K_SYSTEM","R",<Degree(RK),#Coefficients(RK)>,
      "Q",<Degree(QK),#Coefficients(QK)>;
ResK:=CJ!Resultant(RK,QK);
containsMain:=IsDivisibleBy(ResK,mainSR);
print "K_RESULTANT",<TotalDegree(ResK),#Terms(ResK)>,
      "contains_main",containsMain;

KF:=FieldOfFractions(CJ);
KPK<kk>:=PolynomialRing(KF);
qquo,krem:=Quotrem(KPK!QK,KPK!RK);
assert Degree(krem) eq 1;
l0:=Coefficient(krem,0); l1:=Coefficient(krem,1);
krec:=-l0/l1;
print "K_LINEAR_REMAINDER",
      "l0_num",<TotalDegree(Numerator(l0)),#Terms(Numerator(l0))>,
      "l0_den",<TotalDegree(Denominator(l0)),#Terms(Denominator(l0))>,
      "l1_num",<TotalDegree(Numerator(l1)),#Terms(Numerator(l1))>,
      "l1_den",<TotalDegree(Denominator(l1)),#Terms(Denominator(l1))>;
condR:=Evaluate(KPK!RK,krec);
condQ:=Evaluate(KPK!QK,krec);
assert IsDivisibleBy(CJ!Numerator(condR),mainSR);
assert IsDivisibleBy(CJ!Numerator(condQ),mainSR);

coverA:=CJ!(Numerator(l1)*Denominator(l0));
coverB:=CJ!(Numerator(l0)*Denominator(l1));
cg:=GCD(coverA,coverB);
if cg ne 0 and TotalDegree(cg) gt 0 then
    coverA:=ExactQuotient(coverA,cg);
    coverB:=ExactQuotient(coverB,cg);
end if;
print "EXACT_DOUBLE_COVER equation coverA*j^2+coverB=0 on P44";
print " COVER_A",<TotalDegree(coverA),Degree(coverA,sJ),Degree(coverA,rJ),#Terms(coverA)>;
print " COVER_B",<TotalDegree(coverB),Degree(coverB,sJ),Degree(coverB,rJ),#Terms(coverB)>;

// If coverA=coverB=0, the linear k-recovery degenerates.  Audit the
// resulting finite superset directly; emptiness of its Q-variety is enough
// and does not require any denominator saturation.
excG:=GCD(mainSR,GCD(coverA,coverB));
print "EXCEPTIONAL_RECOVERY_GCD",<TotalDegree(excG),#Terms(excG)>;
Iexc:=ideal<CJ|mainSR,coverA,coverB>;
print "EXCEPTIONAL_RECOVERY_BASIS_BEGIN";
time Bexc:=GroebnerBasis(Iexc);
excUnit:=&or[g ne 0 and TotalDegree(g) eq 0:g in Bexc];
print "EXCEPTIONAL_RECOVERY","unit",excUnit,"basis",#Bexc;
if not excUnit then
    dexc:=Dimension(Iexc);
    print "EXCEPTIONAL_RECOVERY_DIMENSION",dexc;
    if dexc eq 0 then
        Clex<se,re>:=PolynomialRing(Q,2,"lex");
        toLex:=hom<CJ->Clex|se,re>;
        Ilex:=ideal<Clex|[toLex(g):g in Bexc]>;
        print "EXCEPTIONAL_RECOVERY_LEX_BEGIN";
        time LexExc:=GroebnerBasis(Ilex);
        print "EXCEPTIONAL_RECOVERY_LEX_SHAPES",
              [<TotalDegree(g),Degree(g,se),Degree(g,re),#Terms(g)>:g in LexExc];
        unis:=[g:g in LexExc|Degree(g,se) eq 0 or Degree(g,re) eq 0];
        print "EXCEPTIONAL_RECOVERY_UNIVARIATES",#unis;
        hasQProjection:=false;
        PRat<tvar>:=PolynomialRing(Q);
        for h in unis do
            if Degree(h,se) eq 0 then
                hu:=PRat!Evaluate(h,[PRat!0,tvar]);
            else
                hu:=PRat!Evaluate(h,[tvar,PRat!0]);
            end if;
            rts:=Roots(hu);
            print " EXCEPTIONAL_UNIVARIATE",<Degree(hu),#Terms(hu)>,"Q_ROOTS",rts;
            if #rts gt 0 then hasQProjection:=true; end if;
        end for;
        print "EXCEPTIONAL_RECOVERY_HAS_Q_PROJECTION",hasQProjection;
    end if;
end if;

print "FIXED_FAST_DONE";
quit;
