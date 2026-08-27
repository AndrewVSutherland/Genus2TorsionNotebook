/* x9join.gp — deep 2-congruence hash-join over the X1(9) family (target
 * [9,9]: ANY 2-congruent pair of 9-torsion curves (3,3)... i.e. (2,2)-glues
 * to a Jacobian with [9]x[9] = [9,9] injected, since the 2-division cubic of
 * a 9-torsion curve is always irreducible (a rational 2-point would give
 * Z/18, dead by Mazur).  Two irreducible cubics define 2-congruent curves
 * iff they cut out isomorphic cubic fields, i.e. same polredabs.
 * Prior art searched only to t-height 20 (lane_x1fam); here we go deep.
 * Usage:  gp -q -s 4G x9join.gp   (set H below or via env)
 * Output: X9COLLISION lines with t-pairs + j-invariants + disc-square flag,
 *         X9CYCLIC lines for square-disc (self-gluable) members.
 */
H = if(getenv("X9H") == 0, 400, eval(getenv("X9H")));
print("X9JOIN height H = ", H);

M = Map();
ncurve = 0;
{
for(qd = 1, H,
  for(pn = -H, H,
    if (gcd(pn, qd) != 1, next);
    my(t = pn/qd);
    if (t == 0 || t == 1, next);
    my(c = t^2*(t-1), b = c*(t^2-t+1));
    if (b == 0, next);
    my(a1 = 1 - c);
    \\ 2-division cubic of [a1,-b,-b,0,0]: 4x^3+b2x^2+2b4x+b6,
    \\ b2 = a1^2-4b, b4 = -a1*b, b6 = b^2; monic integral scaled model:
    my(b2 = a1^2 - 4*b, b4 = -a1*b, b6 = b^2);
    my(g = x^3 + b2*x^2 + 8*b4*x + 16*b6);
    my(m = denominator(content(g)));
    \\ scale roots by m: X^3 + b2*m X^2 + 8 b4 m^2 X + 16 b6 m^3
    my(G = x^3 + b2*m*x^2 + 8*b4*m^2*x + 16*b6*m^3);
    my(K = polredabs(G));
    ncurve++;
    my(jv = ellinit([a1, -b, -b, 0, 0]).j);
    my(dsq = issquare(poldisc(G)));
    if (dsq, print("X9CYCLIC t=", t, " field=", K));
    \\ diamond <2> on X1(9) acts by t -> (t-1)/t (same curve, new point):
    \\ only j-distinct members of a field bucket are genuine congruent pairs
    if (mapisdefined(M, K),
        my(L = mapget(M, K));
        for (i = 1, #L, if (L[i][2] != jv,
            print("X9COLLISION t1=", L[i][1], " t2=", t, " field=", K)));
        mapput(M, K, concat(L, [[t, jv]])),
        mapput(M, K, [[t, jv]]));
  );
  if (qd % 50 == 0, print("PROGRESS qd=", qd, " curves=", ncurve));
);
}
print("X9JOIN_DONE curves=", ncurve, " fields=", #M);
quit;
