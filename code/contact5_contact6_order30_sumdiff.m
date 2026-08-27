//////////////////////////////////////////////////////////////////////
//  Sum/difference form of the simultaneous contact-5/contact-6 cover.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/contact5_contact6_order30_sumdiff.txt";
end if;

Q := Rationals();
R<u,s,q> := PolynomialRing(Q, 3, "grevlex");
K := FieldOfFractions(R);

A := (s+q)/2;
e := (s-q)/2;

F2old := u^6*A - u^3*A^4 - u^6*e + 2*u^3*A^2*e^2 - u^3*e^4
    + 6*u^4*A^2 - 6*u^4*e^2 - 15*u^5 - u^4*A
    + 30*u^3*A^2 + 3*u^4*e - 30*u^3*e^2 - 90*u^4
    + 6*u^2*A^2 - 6*u^2*e^2 - 231*u^3 - u^2*A
    - 3*u^2*e - 90*u^2 - 15*u + A + e;
F3old := -u^3*A^3 + u^3*A^2*e + u^3*A*e^2 - u^3*e^3
    + u^4 + 15*u^3*A + u*A^3 - 15*u^3*e
    + u*A^2*e - u*A*e^2 - u*e^3 + 20*u^3
    - 12*u^2*e - 15*u*A - 15*u*e - 20*u - 1;

F2 := R!Numerator(K!F2old);
F3 := R!Numerator(K!F3old);

out := Open(output_file, "w");
fprintf out, "# Sum/difference variables: s=A+e, q=A-e\n\n";
fprintf out, "F2sd = %o\n", F2;
fprintf out, "F3sd = %o\n", F3;
fprintf out, "F2 factors:\n";
for pair in Factorization(F2) do
    fprintf out, "  exponent %o: %o\n", pair[2], pair[1];
end for;
fprintf out, "F3 factors:\n";
for pair in Factorization(F3) do
    fprintf out, "  exponent %o: %o\n", pair[2], pair[1];
end for;
fprintf out, "gcd = %o\n", GCD(F2,F3);

res_q := Resultant(F2, F3, q);
fprintf out, "\nResultant eliminating q: degree %o total_degree %o terms %o\n",
        Degree(res_q), TotalDegree(res_q), #Terms(res_q);
fprintf out, "factorization:\n";
for pair in Factorization(res_q) do
    fprintf out, "  exponent %o degree %o total_degree %o terms %o: %o\n",
            pair[2], Degree(pair[1]), TotalDegree(pair[1]), #Terms(pair[1]), pair[1];
end for;

res_s := Resultant(F2, F3, s);
fprintf out, "\nResultant eliminating s: degree %o total_degree %o terms %o\n",
        Degree(res_s), TotalDegree(res_s), #Terms(res_s);
fprintf out, "factorization:\n";
for pair in Factorization(res_s) do
    fprintf out, "  exponent %o degree %o total_degree %o terms %o: %o\n",
            pair[2], Degree(pair[1]), TotalDegree(pair[1]), #Terms(pair[1]), pair[1];
end for;

delete out;
print "Wrote", output_file;
quit;
