//////////////////////////////////////////////////////////////////////
// Exact verification of general-5 survivors on the compact M(12) chart.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=8; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q:=Rationals(); Z:=Integers(); P<x>:=PolynomialRing(Q); PT<T>:=PolynomialRing(Q);

if assigned BNum and assigned BDen and assigned WNum and assigned WDen then
    if Type(BNum) eq MonStgElt then BNum:=StringToInteger(BNum); end if;
    if Type(BDen) eq MonStgElt then BDen:=StringToInteger(BDen); end if;
    if Type(WNum) eq MonStgElt then WNum:=StringToInteger(WNum); end if;
    if Type(WDen) eq MonStgElt then WDen:=StringToInteger(WDen); end if;
    candidates := [<Q!BNum/Q!BDen,Q!WNum/Q!WDen>];
else
    // The ten survivors of the full height-200 scan through p=79.
    candidates := [
      <Q!-61/23,Q!-87/187>, <Q!-154/43,Q!-52/97>,
      <Q!183/67,Q!-51/86>, <Q!141/70,Q!187/150>,
      <Q!79/74,Q!117/145>, <Q!91/114,Q!-59/142>,
      <Q!115/134,Q!-160/143>, <Q!124/177,Q!-143/36>,
      <Q!-32/183,Q!-76/141>, <Q!-193/198,Q!-26/125>
    ];
end if;

function CompactQuintic(b,w)
    L:=b+(2*b-1)*x;
    A:=x+w*(1+b*x);
    return P!(L*(L*A^2+4*b*(1+x)^2*(w*L-x^2)));
end function;

function IntegralModel(f)
    d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return P!(d^2*f),d;
end function;

function FrobeniusPolynomial(C,p)
    ef:=EulerFactor(C,p); d:=Degree(ef);
    return &+[Q!Coefficient(ef,i)*T^(d-i):i in [0..d]];
end function;

function FullSimplicityCertificate(C, bound)
    for p in [ell:ell in PrimesUpTo(bound)|ell ge 3] do
        try
            Phi:=FrobeniusPolynomial(C,p);
            fac:=Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then continue; end if;
            Gal:=GaloisGroup(Phi); desc:="";
            try desc:=TransitiveGroupDescription(Gal); catch e desc:="unknown"; end try;
            if Order(Gal) eq 8 and desc eq "D(4)" then
                return true,"D4",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    for p in [ell:ell in PrimesUpTo(bound)|ell ge 3] do
        try
            Phi:=FrobeniusPolynomial(C,p);
            fac:=Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then continue; end if;
            K<pi>:=NumberField(Phi); ok:=true;
            for n in [2..12] do
                if Degree(MinimalPolynomial(pi^n)) lt 4 then ok:=false; break; end if;
            end for;
            if ok then return true,"root_power",p,Phi; end if;
        catch e
            continue;
        end try;
    end for;
    return false,"none",0,PT!0;
end function;

seen:={};
for bw in candidates do
    b:=bw[1]; w:=bw[2]; f:=CompactQuintic(b,w);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        print " SINGULAR_OR_DEGREE_BOUNDARY";
        continue;
    end if;
    fI,d:=IntegralModel(f);
    key:=Sprint(fI);
    print "BW_CANDIDATE","b",b,"w",w;
    print " f",f;
    print " f_integral",fI,"yscale",d;
    C:=HyperellipticCurve(fI); J:=Jacobian(C);
    eta:=(1-b)*(w*(1-b)-1);
    assert Evaluate(f,-1) eq eta^2;
    D12:=J![x+1,d*eta];
    print " D12_order",Order(D12),"D12",D12;
    if key in seen then print " DUPLICATE_INTEGRAL_MODEL"; continue; end if;
    Include(~seen,key);
    G,phi:=TorsionSubgroup(J); inv:=Invariants(G);
    print " TORSION",inv,"order",#G;
    for i in [1..Ngens(G)] do
        gi:=G.i; Pi:=phi(gi);
        print " GENERATOR",i,"order",Order(gi),"mumford",Pi;
        if Order(gi) mod 60 eq 0 then
            P60:=(Order(gi) div 60)*Pi;
            print " P60",P60,"order",Order(P60);
            print " P5",12*P60,"order",Order(12*P60);
            print " P12",5*P60,"order",Order(5*P60);
        end if;
    end for;
    simple,stype,sp,Phi:=FullSimplicityCertificate(C,127);
    print " GEOMETRICALLY_SIMPLE",simple,"type",stype,"prime",sp;
    if simple then print " FROBENIUS",Phi; end if;
    print " CYCLIC60_CERTIFIED",inv eq [60] and simple;
end for;
print "GENERAL5_VERIFY_DONE";
quit;
