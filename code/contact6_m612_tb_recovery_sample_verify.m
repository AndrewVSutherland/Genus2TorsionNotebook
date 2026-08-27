//////////////////////////////////////////////////////////////////////
// Exact sample verification on the rational recovery-boundary family.
// The parameter u=2 gives s=19/30 and the two rational a-branches below.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); P<x>:=PolynomialRing(Q);
PF<T>:=PolynomialRing(Q);

s:=Q!19/30;
r:=Q!608/225;
K:=Q!153664/50625;
assert K eq (Q!392/225)^2;
b:=2*s^2-3;
alist:=[Q!1601/2025,Q!48791/1875];

function IntegralModel(f)
    d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return P!(d^2*f);
end function;

function FrobeniusPolynomial(C,p)
    ef:=EulerFactor(C,p); d:=Degree(ef);
    return &+[Q!Coefficient(ef,i)*T^(d-i):i in [0..d]];
end function;

function FullSimplicityCertificate(C)
    plist:=[p:p in PrimesUpTo(151)|p ge 3];
    for p in plist do
        try
            Phi:=FrobeniusPolynomial(C,p);
            fac:=Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            Gal:=GaloisGroup(Phi);
            desc:="unknown";
            try desc:=TransitiveGroupDescription(Gal);
            catch e2 desc:="unknown"; end try;
            if Order(Gal) eq 8 and desc eq "D(4)" then
                return true,"D4",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    for p in plist do
        try
            Phi:=FrobeniusPolynomial(C,p);
            fac:=Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            Kf<pi>:=NumberField(Phi);
            if &and[Degree(MinimalPolynomial(pi^n)) eq 4:n in [2..12]] then
                return true,"root_power",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,"none",0,PF!0;
end function;

for a in alist do
    h:=1+a*x+b*x^2+x^3;
    f:=h^2-(x-1)^6;
    print "SAMPLE","s",s,"r",r,"K",K,"a",a,"b",b;
    print "SMOOTH",Degree(f),Discriminant(f) ne 0,
          "FACTORS",Factorization(f);
    fI:=IntegralModel(f);
    C:=HyperellipticCurve(fI); J:=Jacobian(C);
    TG,mp:=TorsionSubgroup(J);
    print "TORSION",Invariants(TG);
    simple,method,pcert,Phi:=FullSimplicityCertificate(C);
    print "FULL_SIMPLE",simple,method,pcert,Phi;
    try print "AUT_ORDER",#AutomorphismGroup(C);
    catch e print "AUT_FAILED"; end try;
    for p in [5,7,11,13,17,19,23,29,31,37,41,43] do
        try
            fp:=ChangeRing(fI,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            Lp:=LPolynomial(ChangeRing(C,GF(p)));
            print "LP",p,Lp,Factorization(Lp);
        catch e
            print "LP_FAILED",p;
        end try;
    end for;
end for;
quit;
