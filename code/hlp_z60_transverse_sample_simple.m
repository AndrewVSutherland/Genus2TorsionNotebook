//////////////////////////////////////////////////////////////////////
// A geometric-simplicity certificate for the sample t=1 on the
// transverse line F_t=F_0+t*x.  This sample is not asserted to retain
// rational 60-torsion; it proves that the selected line actually reaches
// the geometrically simple locus.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);
Q:=Rationals(); P<x>:=PolynomialRing(Q); PT<T>:=PolynomialRing(Q);

f:=-46250000*x^6+1761500625*x^4-22332312000*x^2+x+94277468160;
C:=HyperellipticCurve(f);

function FrobeniusPolynomial(C,p)
    ef:=EulerFactor(C,p); d:=Degree(ef);
    return &+[Q!Coefficient(ef,i)*T^(d-i):i in [0..d]];
end function;

print "HLP_Z60_TRANSVERSE_SAMPLE_SIMPLE";
print "t",1,"f",f;
for p in [ell:ell in PrimesUpTo(97)|ell ge 3] do
    try
        Phi:=FrobeniusPolynomial(C,p);
        fac:=Factorization(Phi);
        if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
            continue;
        end if;
        K<pi>:=NumberField(Phi); ok:=true;
        degrees:=[];
        for n in [1..12] do
            dn:=Degree(MinimalPolynomial(pi^n));
            Append(~degrees,dn);
            if dn lt 4 then ok:=false; end if;
        end for;
        if ok then
            print "ABSOLUTELY_SIMPLE_REDUCTION_CERTIFICATE prime",p;
            print "Frobenius_polynomial",Phi;
            print "power_degrees_1_to_12",degrees;
            print "CONCLUSION geometrically_simple_over_Q true";
            quit;
        end if;
    catch e
        continue;
    end try;
end for;
print "NO_CERTIFICATE_FOUND_THROUGH_97";
quit;
