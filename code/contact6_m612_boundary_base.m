//////////////////////////////////////////////////////////////////////
// Exact base boundary geometry for the contact-6 [6,6] core.
//
// The normalized family is
//   f=x*B*C,
//   B=(b+3)x^2+(a-3)x+2,
//   C=2x^2+(b-3)x+(a+3).
//
// This script simplifies the advertised RR boundary, homogenizes every
// boundary support in [A:B:T]=[a:b:1], and enumerates its complete
// intersection stratification over F_5 and F_7.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();
R<a,b>:=PolynomialRing(Q,2);
P<x>:=PolynomialRing(R);

BF:=(b+3)*x^2+(a-3)*x+2;
CF:=2*x^2+(b-3)*x+(a+3);
f:=x*BF*CF;
DB:=(a-3)^2-8*(b+3);
DC:=(b-3)^2-8*(a+3);
DeltaR:=(a+3)*(b+3)-4;
RR:=DeltaR^2-(b^2-2*a-3)*(a^2-2*b-3);

assert Resultant(BF,CF) eq 2*(a+b+2)^3;
assert RR eq 2*(a+b+2)^3;
disc:=Discriminant(f);
disc_support:=(a+3)^2*(a+b+2)^6*DB*DC;
assert IsDivisibleBy(disc,disc_support);
disc_unit:=ExactQuotient(disc,disc_support);
assert Degree(disc_unit) eq 0;

print "CONTACT6_M612_BOUNDARY_BASE";
print "RESULTANT_BC",Resultant(BF,CF);
print "RR_SIMPLIFIED",RR;
print "DISC_FACTOR_SUPPORT",
      "(a+3)^2*(a+b+2)^6*DB*DC";
print "DISC_UNIT",disc_unit;
print "NOTE b+3 is the degree-loss divisor and is not detected by the";
print "     specialized degree-5 discriminant polynomial.";

// Projective coordinates [A:B:T].  The five affine supports and the
// parameter-infinity line are all defined integrally at p=5,7.
H<A,B,T>:=PolynomialRing(Integers(),3);
Lb:=B+3*T;
La:=A+3*T;
LR:=A+B+2*T;
DBh:=(A-3*T)^2-8*T*(B+3*T);
DCh:=(B-3*T)^2-8*T*(A+3*T);
Linf:=T;
boundaries:=[Lb,La,LR,DBh,DCh,Linf];
labels:=["b+3","a+3","a+b+2 (= RR support)","DB","DC","infinity"];

print "PROJECTIVE_SUPPORTS";
for i in [1..#labels] do print labels[i],boundaries[i]; end for;
print "RATIONAL_PARAMETRIZATIONS [u:v] -> [A:B:T]";
print " b+3: [u:-3*v:v]";
print " a+3: [-3*v:u:v]";
print " a+b+2: [u:-u-2*v:v]";
print " DB: [8*u*v:u^2-6*u*v-15*v^2:8*v^2]";
print " DC: [u^2-6*u*v-15*v^2:8*u*v:8*v^2]";
print " infinity: [u:v:0]";
print "INFINITY_SPECIAL_POINTS",
      "[1:0:0] on b+3 and tangent DC;",
      "[0:1:0] on a+3 and tangent DB;",
      "[1:-1:0] on a+b+2.";

function MaskLabels(mask,labels)
    ans:=[];
    bit:=1;
    for i in [1..#labels] do
        if (mask div bit) mod 2 eq 1 then Append(~ans,labels[i]); end if;
        bit*:=2;
    end for;
    return ans;
end function;

for p in [5,7] do
    k:=GF(p);
    P2:=ProjectiveSpace(k,2);
    CR:=CoordinateRing(P2);
    bmods:=[Evaluate(ChangeRing(g,k),[CR.1,CR.2,CR.3])
            :g in boundaries];
    strata:=AssociativeArray();
    points_by_mask:=AssociativeArray();
    union_count:=0;
    for pt in Points(P2) do
        coords:=Eltseq(pt);
        mask:=0; bit:=1;
        for g in bmods do
            if Evaluate(g,coords) eq 0 then mask+:=bit; end if;
            bit*:=2;
        end for;
        if mask eq 0 then continue; end if;
        union_count+:=1;
        if IsDefined(strata,mask) then
            strata[mask]+:=1;
            Append(~points_by_mask[mask],
                   <Integers()!coords[1],Integers()!coords[2],
                    Integers()!coords[3]>);
        else
            strata[mask]:=1;
            points_by_mask[mask]:=[
                <Integers()!coords[1],Integers()!coords[2],
                 Integers()!coords[3]>];
        end if;
    end for;
    print "PRIME",p,"P2_POINTS",#Points(P2),"BOUNDARY_UNION",union_count;
    for mask in Sort(Setseq(Keys(strata))) do
        labs:=MaskLabels(mask,labels);
        print " STRATUM",mask,labs,"COUNT",strata[mask];
        if #labs ge 2 then print "  POINTS",points_by_mask[mask]; end if;
    end for;
end for;

print "CONTACT6_M612_BOUNDARY_BASE_DONE";
quit;
