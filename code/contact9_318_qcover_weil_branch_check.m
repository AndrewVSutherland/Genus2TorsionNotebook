// Assign the two independent t=4 mod 7 order-3 branches to a Weil orbit.

F := GF(7);
P<x> := PolynomialRing(F);

// Marked line 3*D9.  Here G=H/h3 and s=1/h3, so z=s^2.
q0 := x^2 + x + 4;
G0 := x^3 + 6*x^2 + 6*x + 2;
s0 := F!3;
assert s0^2 eq 2;

candidates := [
    <x^2 + 2*x + 3, x^3 + 4*x^2 + 5*x + 3, F!3, "branch_(2,3,2)">,
    <x^2 + 3*x + 3, x^3 + 5*x^2 + 5*x + 2, F!1, "branch_(3,3,1)">
];

print "contact9 [3,18] t=4 mod 7 Weil-orbit check";
print "marked", q0, G0, "s", s0;
for item in candidates do
    q := item[1];
    G := item[2];
    s := item[3];
    label := item[4];
    lambda := s/s0;
    value_on_marked := Resultant(q0, G - lambda*G0);
    value_on_candidate_scaled := Resultant(q, G - lambda*G0);
    pairing_equation := lambda^2*value_on_marked - value_on_candidate_scaled;
    print label,
          "lambda", lambda,
          "marked_resultant", value_on_marked,
          "candidate_resultant", value_on_candidate_scaled,
          "pairing_equation", pairing_equation,
          "orthogonal", pairing_equation eq 0;
    assert pairing_equation eq 0;
end for;

print "DONE both independent branches lie on the degree-12 orthogonal orbit";
quit;
