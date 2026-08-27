//////////////////////////////////////////////////////////////////////
// Direct point-count profiles of the seven pair-scaling branch families.
//////////////////////////////////////////////////////////////////////
SetColumns(0); Z:=Integers();
families:=[
 <"P1",[1,55,99,125]>,<"P2",[20,225,304,380]>,
 <"P3",[13,276,624,828]>,<"P4",[41,256,800,1312]>,
 <"P5",[2,46,292,1679]>,<"P6",[18,158,711,1764]>,
 <"P7",[50,791,2800,5650]>
];
primes:=[11,13,17,19,23,29,31];
for fam in families do
    name:=fam[1]; z:=fam[2]; a:=z[1];b:=z[2];c:=z[3];d:=z[4];
    for ell in primes do
        F:=GF(ell); R<x>:=PolynomialRing(F);
        boundary:=0; good:=0; target:=0; residues:=[];
        for li in [0..ell] do
            if li eq ell then vals:=[F!a,F!0,F!0,F!d]; label:="inf";
            else lam:=F!li; vals:=[lam*a,F!b,F!c,lam*d]; label:=IntegerToString(li);
            end if;
            f:=x*&*[x+q^2:q in vals];
            if Degree(f) ne 5 or Discriminant(f) eq 0 then boundary+:=1; continue; end if;
            good+:=1; n:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
            if n mod 3 eq 0 then target+:=1; Append(~residues,label); end if;
        end for;
        print "PAIRFIBER_PROFILE",name,"p",ell,"boundary",boundary,
              "good",good,"target",target,"residues",residues;
    end for;
end for;
quit;
