//////////////////////////////////////////////////////////////////////
// Exact arithmetic audit of the rational point u=12 on the P8 R3
// component: e=-200/409, hence (a,b)=(-409/200,0).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
run_known_seed:=false; NoMain:=true;
load "code/contact6_m612_dual_class_exact.m";

Q:=Rationals(); P<x>:=PolynomialRing(Q);

function IntegralSquareScaleLocal(f)
    d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return P!(d^2*f);
end function;

function SimpleCertificateLocal(f)
    C:=HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67] do
        fp:=ChangeRing(f,GF(p));
        if Degree(fp) ne Degree(f) or Discriminant(fp) eq 0 then continue; end if;
        Lp:=LPolynomial(ChangeRing(C,GF(p)));
        fac:=Factorization(Lp);
        if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
            return true,p,Lp;
        end if;
    end for;
    return false,0,PolynomialRing(Integers())!0;
end function;

e0:=Q!-200/409; a:=1/e0; b:=Q!0;
assert a eq Q!-409/200;
B:=(b+3)*x^2+(a-3)*x+2;
C:=2*x^2+(b-3)*x+(a+3);
f:=x*B*C;
assert Degree(f) eq 5 and Discriminant(f) ne 0;
fI:=IntegralSquareScaleLocal(f);
Js:=Jacobian(HyperellipticCurve(fI));
Gs,mps:=TorsionSubgroup(Js); source_inv:=Invariants(Gs);
source_simple,source_p,source_L:=SimpleCertificateLocal(fI);

Delta,R1,R2,R3,DB,DC:=DualData(a,b);
g:=Delta*R1*R2*R3;
assert Degree(g) eq 6 and Discriminant(g) ne 0;
gI:=IntegralSquareScaleLocal(g);
Jd:=Jacobian(HyperellipticCurve(gI));
Gd,mpd:=TorsionSubgroup(Jd); dual_inv:=Invariants(Gd);
dual_simple,dual_p,dual_L:=SimpleCertificateLocal(gI);
classes:=[Jd![q/LeadingCoefficient(q),0]:q in [R1,R2,R3]];
halves:=[IsDivisibleBy(T,2):T in classes];

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P8_AUDIT";
print "POINT","e",e0,"a",a,"b",b;
print "SOURCE_SQUAREFREE",IsSquarefree(f),"TORSION",source_inv,
      "TWO_3_DIRECTIONS",#[n:n in source_inv|n mod 3 eq 0] ge 2;
print "SOURCE_SIMPLE_CERTIFICATE",source_simple,"PRIME",source_p,"LPOLY",source_L;
print "DUAL_SQUAREFREE",IsSquarefree(g),"TORSION",dual_inv,
      "TWO_3_DIRECTIONS",#[n:n in dual_inv|n mod 3 eq 0] ge 2;
print "DUAL_HALVES_R1_R2_R3",halves;
print "DUAL_SIMPLE_CERTIFICATE",dual_simple,"PRIME",dual_p,"LPOLY",dual_L;
assert source_inv eq [2,2,6];
assert dual_inv eq [2,12];
assert halves eq [false,false,true];
assert #[n:n in source_inv|n mod 3 eq 0] eq 1;
assert #[n:n in dual_inv|n mod 3 eq 0] eq 1;
print "VERDICT","R3 half is genuine, but neither source nor dual has rational [3,3]";
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P8_AUDIT_DONE";
quit;
