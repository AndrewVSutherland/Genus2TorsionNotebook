//////////////////////////////////////////////////////////////////////
// Independent verification of the new [2,2,2,12] order-96 record.
//
// Checks:
//   1. rebuild the collaborator's generalized hyperelliptic model;
//   2. compare its simplified model with the reported sextic;
//   3. compute the exact rational torsion subgroup;
//   4. produce several good-prime root-power Frobenius certificates
//      for geometric simplicity.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

f0 := [
    830742747091037849,
    32014154874551031,
    -530648977741620,
    -15854483576121,
    150572203590,
    737595570,
    756900
];
h0 := [1,0,1];
C0 := HyperellipticCurve([f0,h0]);
Cs := SimplifiedModel(C0);
fs, hs := HyperellipticPolynomials(Cs);

f_reported :=
    3027600*x^6 + 2950382280*x^5 + 602288814361*x^4
    - 63417934304484*x^3 - 2122595910966478*x^2
    + 128056619498204124*x + 3322970988364151397;

print "INPUT_CURVE", C0;
print "SIMPLIFIED_CURVE", Cs;
print "SIMPLIFIED_POLYNOMIALS", fs, hs;
print "REPORTED_SEXTIC_MATCH", hs eq 0 and fs eq f_reported;
assert hs eq 0 and fs eq f_reported;
assert Discriminant(fs) ne 0;

J := Jacobian(Cs);
T, phi := TorsionSubgroup(J);
invs := Invariants(T);
print "TORSION_INVARIANTS", invs;
print "TORSION_ORDER", #T;
assert invs eq [2,2,2,12];
assert #T eq 96;

function CountCurve(fp)
    Fq := BaseRing(Parent(fp));
    count := 0;
    for xx in Fq do
        value := Evaluate(fp,xx);
        if value eq 0 then
            count +:= 1;
        elif IsSquare(value) then
            count +:= 2;
        end if;
    end for;
    if IsSquare(LeadingCoefficient(fp)) then count +:= 2; end if;
    return count;
end function;

function RootPowerWitness(f,pp)
    if Z!LeadingCoefficient(f) mod pp eq 0 then
        return false, Parent(f)!0, [];
    end if;
    if Z!Numerator(Discriminant(f)) mod pp eq 0 then
        return false, Parent(f)!0, [];
    end if;
    Fp := GF(pp); Pp := PolynomialRing(Fp);
    fp := Pp![Fp!c : c in Coefficients(f)];
    if Degree(fp) lt 5 or not IsSquarefree(fp) then
        return false, Parent(f)!0, [];
    end if;
    Fp2 := GF(pp^2); Pp2 := PolynomialRing(Fp2);
    fp2 := Pp2![Fp2!c : c in Coefficients(f)];
    a1 := pp + 1 - CountCurve(fp);
    a2 := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;
    R<TT> := PolynomialRing(Q);
    chi := TT^4 - a1*TT^3 + a2*TT^2 - a1*pp*TT + pp^2;
    Lp := LPolynomial(HyperellipticCurve(fp));
    assert [Z!Coefficient(Lp,i) : i in [0..4]]
           eq [1,-a1,a2,-a1*pp,pp^2];
    if not IsIrreducible(chi) then return false, chi, []; end if;
    K<pi> := NumberField(chi);
    degrees := [Degree(MinimalPolynomial(pi^n)) : n in [2..12]];
    return &and[d eq 4 : d in degrees], chi, degrees;
end function;

// Independent prime-to-bad-reduction bound for exactness of the torsion
// order: rational torsion injects into both finite Jacobians.
ok31, chi31, deg31 := RootPowerWitness(f_reported,31);
ok37, chi37, deg37 := RootPowerWitness(f_reported,37);
order31 := Z!Evaluate(chi31,1);
order37 := Z!Evaluate(chi37,1);
print "REDUCTION_TORSION_BOUND", "p31_order", order31,
      "p37_order", order37, "gcd", GCD(order31,order37);
assert order31 eq 864 and order37 eq 1248;
assert GCD(order31,order37) eq 96;

witnesses := [];
for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
    ok, chi, degrees := RootPowerWitness(f_reported,pp);
    if chi ne 0 then
        printf "FROBENIUS p=%o chi=%o degrees_n2_to_12=%o witness=%o\n",
            pp, chi, degrees, ok;
    end if;
    if ok then Append(~witnesses,<pp,chi>); end if;
end for;

print "GEOMETRIC_SIMPLICITY_WITNESSES", witnesses;
assert #witnesses ge 2;
print "RECORD_22212_ORDER96_VERIFIED";
quit;
