//////////////////////////////////////////////////////////////////////
// Validate an explicit Weil-pairing test for the marked 3-line versus
// a cubic-contact 3-class.  The test is later used on p=7 weighted
// initial branches to distinguish the orthogonal degree-12 support.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
k:=GF(31); P<x>:=PolynomialRing(k);
tau:=k!7;
t:=4*(tau^2+tau-6)/(tau^2+6);
e:=-(k!25/k!3)*t^2/(t^4-k!25*t^2+k!1250/k!3);
a:=1/e;
h6:=1+a*x+x^3; f:=h6^2-(x-1)^6;
assert Degree(f) eq 5 and Discriminant(f) ne 0;
C:=HyperellipticCurve(f); J:=Jacobian(C);
uP:=(x-1)^2;
P3:=J![uP,h6 mod uP];
assert Order(P3) eq 3;

checked:=0; orthogonal:=0; nonorthogonal:=0; exact_equal:=0;
inverse_equal:=0; pair_values:={k|};
for L in k do
    if L eq 0 then continue; end if;
    M:=L^2;
    for U in k do for V in k do
        if V eq 0 or U^2-4*V^2 eq 0 then continue; end if;
        N:=(3*U+6*M)/2;
        R:=(3*U^2+3*V^2+(2/e-15)*M-N^2)/2;
        F3:=2*V^3+2*N*R-U^3-6*U*V^2-22*M;
        F2:=R^2+2*N*V^3-3*U^2*V^2-3*V^4-(1/e^2-15)*M;
        F1:=2*R*V^3-3*U*V^4-(2/e+6)*M;
        if F1 ne 0 or F2 ne 0 or F3 ne 0 then continue; end if;
        q:=x^2+U*x+V^2;
        H:=x^3+N*x^2+R*x+V^3;
        if Degree(GCD(q,f)) gt 0 then continue; end if;
        Q3:=J![q,(H/L) mod q];
        assert Order(Q3) eq 3;
        wp:=WeilPairing(P3,Q3,3);
        num:=Resultant(q,H-L*h6);
        den:=L^2*Resultant(uP,L*h6-H);
        if den eq 0 then continue; end if;
        explicit:=num/den;
        assert (wp eq 1) eq (explicit eq 1);
        if wp eq explicit then exact_equal+:=1; end if;
        if wp eq explicit^-1 then inverse_equal+:=1; end if;
        Include(~pair_values,wp); checked+:=1;
        if wp eq 1 then orthogonal+:=1; else nonorthogonal+:=1; end if;
    end for; end for;
end for;

print "CONTACT6_M612_P8_WEIL_PAIRING";
print "PRIME",31,"TAU",7,"E",e,"CHECKED",checked,
      "ORTHOGONAL",orthogonal,"NONORTHOGONAL",nonorthogonal,
      "EXACT_EQUAL",exact_equal,"INVERSE_EQUAL",inverse_equal,
      "PAIR_VALUES",pair_values;
assert checked gt 0 and orthogonal gt 0 and nonorthogonal gt 0;
assert exact_equal eq checked or inverse_equal eq checked;
print "FORMULA",
      "Res(q,H-L*h6)/(L^2*Res((x-1)^2,L*h6-H)); orthogonal iff 1";
print "CONTACT6_M612_P8_WEIL_PAIRING_DONE";
quit;
