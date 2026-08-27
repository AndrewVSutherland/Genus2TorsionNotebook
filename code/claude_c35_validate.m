// Validation for C35 test run (target Z/35, claude_top10_10_35)
// (1) Elkies LuCaNT A_1(5) universal 5-torsion identity at sample (q0,q1,q2)
// (2) repo contact-7 family [28] curve: torsion + simplicity certificate
QQ := Rationals();
R<x> := PolynomialRing(QQ);

print "=== Validation 1: Elkies A_1(5) universal 5-torsion ===";
for triple in [[1,1,1],[2,-1,3],[1/2,3,-2],[-3,2/5,7]] do
    q0 := QQ!triple[1]; q1 := QQ!triple[2]; q2 := QQ!triple[3];
    Q  := q2*x^2 + q1*x + q0;
    Qp := Q - x;              // Q' = Q - L*L' with L=x, L'=1
    H  := Qp - x*Q;           // L'Q' - LQ
    f  := Q^2*Qp;
    fz := H^2 + 4*f;          // z^2 = fz model
    if Degree(fz) lt 5 or not IsSquarefree(fz) then
        print triple, "DEGENERATE";
        continue;
    end if;
    C := HyperellipticCurve(f, H);
    J := Jacobian(C);
    a := Q/LeadingCoefficient(Q);
    T := J![a, R!0];
    print triple, "genus", Genus(C), "order_T", Order(T);
end for;

print "=== Validation 2: contact-7 family [28] curve ===";
f28 := 4*x^5 + 21*x^4 - 70*x^3 + 79*x^2 - 42*x + 9;
C28 := HyperellipticCurve(f28);
J28 := Jacobian(C28);
T28 := TorsionSubgroup(J28);
print "torsion invariants:", Invariants(T28);
// marked 7-torsion class J![x-1, h(1)] with a=11/2,b=-7/2 scaling: repo says J![x-1,-1]
D7 := J28![x-1, R!-1];
print "marked class order:", Order(D7);
// simplicity certificate at p=5
C5 := ChangeRing(C28, GF(5));
Lp := LPolynomial(C5);
print "L_5 =", Lp;
cp := R!Reverse(Coefficients(Lp));  // charpoly of Frobenius (monic quartic)
print "charpoly irreducible:", IsIrreducible(cp);
Rxy<X,Y> := PolynomialRing(QQ,2);
res := Resultant(Evaluate(cp,X), Y - X^12, X);
m12 := UnivariatePolynomial(res);
print "12th-power transform degree:", Degree(m12), "irreducible:", IsIrreducible(m12);

print "=== Validation 3: Howe order-70 split curve must be bad at 3 ===";
f70 := 22*x^5 + (697/144)*x^4 - (645/4)*x^3 + (1045/4)*x^2 - 162*x + 36;
// integral model: y -> y/144, x same: 144^2*f70
g70 := R!(144^2*f70);
print "disc valuation at 3 of integral model:", Valuation(Discriminant(g70),3);
print "(Weil: good reduction at 3 impossible since 35 | #J(F_3) needs #J=35<=55; 70>55)";
quit;
