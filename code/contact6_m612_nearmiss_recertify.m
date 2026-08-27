// Re-certify the [2,12] near-miss (a,b)=(-409/200,0) dual curve with a
// VALID (non-even chi) D4/root-power certificate.  The contact-6 dual:
// rebuild from the P8 data: e=-200/409 fiber, a=1/e, b=0 on contact-6:
// h6 = 1 + a*x + b*x^2 + x^3, f = h6^2 - (x-1)^6 with (a,b)=(-409/200,0).
SetColumns(0);
Q := Rationals(); Z := Integers(); P<x> := PolynomialRing(Q);
a := -409/200; b := 0;
h6 := 1 + a*x + b*x^2 + x^3;
f := h6^2 - (x-1)^6;
L := LCM([Denominator(c) : c in Coefficients(f)]);
fInt := P!(L^2*f);
printf "fInt = %o\n", fInt;
J := Jacobian(HyperellipticCurve(fInt));
inv := Invariants(TorsionSubgroup(J));
printf "torsion = %o\n", inv;
// simplicity: root-power certificate, REJECTING even chi
RT := PolynomialRing(Q); T := RT.1;
dsc := Discriminant(fInt);
function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
    if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if;
    return cnt;
end function;
found := false;
for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73] do
    if (Z!LeadingCoefficient(fInt)) mod pp eq 0 then continue; end if;
    if (Z!Numerator(dsc)) mod pp eq 0 then continue; end if;
    PF := PolynomialRing(GF(pp));
    fp := PF![GF(pp)!co : co in Coefficients(fInt)];
    if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
    PF2 := PolynomialRing(GF(pp^2));
    fp2 := PF2![GF(pp^2)!co : co in Coefficients(fInt)];
    a1 := pp + 1 - CountCurve(fp);
    a2 := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;
    chi := T^4 - a1*T^3 + a2*T^2 - a1*pp*T + pp^2;
    if a1 eq 0 then continue; end if;   // even chi cannot certify: skip
    if not IsIrreducible(chi) then continue; end if;
    K := NumberField(chi); pi := K.1; drop := false;
    for nn in [2..12] do
        if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
    end for;
    if not drop then
        printf "VALID CERTIFICATE at p=%o : chi = %o (a1=%o, odd part nonzero)\n", pp, chi, a1;
        found := true; break;
    end if;
end for;
if not found then print "NO valid odd certificate through p=73 -- simplicity UNRESOLVED"; end if;
print "DONE";
quit;
