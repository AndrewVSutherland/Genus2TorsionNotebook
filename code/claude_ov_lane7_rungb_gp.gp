\\ claude_ov_lane7_rungb_gp.gp -- Lane 7 (overnight 2026-07-25).
\\
\\ INDEPENDENT (PARI) cross-check of RUNG (b) from code/claude_ov_lane7_rungs.m,
\\ plus the FIRST FAILING PRIME that the Magma run did not record.
\\
\\ Rung (b) targets order 168 = [2,2,42]: an INDEPENDENT rational 3-torsion
\\ class on a contact-7 three-root member.  Necessary condition: 3 | #J(F_p) at
\\ EVERY good prime p != 3.  Here #J(F_p) = chi_p(1) with chi_p the Frobenius
\\ characteristic polynomial from hyperellcharpoly -- no Magma, no Jacobian
\\ group computation.
\\
\\ CONTROL: the [6,6] witness has full rational 3-torsion, so it must pass at
\\ 100% of good primes.
\\
\\ Markers: RUNGBGP / CONTROLGP / RUNGBGPDONE

G(v) = -(v^5 - v^3 - v^2/2)/(v+1)^2;

trips = [ [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5], \
          [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17], \
          [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289], \
          [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], \
          [-1/2, -13/49, 13/50] ];

PMAX = 500;

test3(F, lbl) = {
  my(D = poldisc(F), lc = pollead(F), tot = 0, pass = 0, first = 0, skipped = 0, nJ);
  forprime(p = 11, PMAX,
    if (lc % p == 0, next);                 \\ leading coeff vanishes: bad model at p
    if (numerator(D) % p == 0, next);       \\ bad reduction
    \\ everything from here is wrapped: hyperellcharpoly raises on any model it
    \\ cannot handle mod p, and the raise escapes a narrower iferr
    nJ = iferr(
          my(Fp = Mod(1,p) * F);
          if (poldegree(Fp) != 5 || !issquarefree(Fp), 0,
              subst(hyperellcharpoly(Fp), x, 1)),
         E, skipped++; 0);
    if (nJ == 0, next);
    tot++;
    if (nJ % 3 == 0, pass++, if (first == 0, first = p));
  );
  printf("%s : 3 | #J(F_p) at %d/%d good primes (%d%%)  first failing prime = %d  [skipped %d]\n",
         lbl, pass, tot, if(tot, 100*pass\tot, 0), first, skipped);
};

{
print("---- CONTROL: the [6,6] witness (full rational 3-torsion) ----");
test3(11389248*x^5 - 18252000*x^4 + 42399396*x^3 - 10288044*x^2 + 29659500*x, "CONTROLGP [6,6]");

print("---- the eleven contact-7 three-root points ----");
for(i = 1, #trips,
  my(s = trips[i][1], t = trips[i][2], u = trips[i][3]);
  my(c4 = (G(s)-G(t))/(s^2-t^2));
  my(c0 = G(s) - c4*s^2);
  my(b = c4 - 2, a = 9/2 - c0 - c4);
  my(h = 1 - 7/2*x + a*x^2 + b*x^3);
  my(f = (h^2 + (x-1)^7) \ x^2);
  my(d = lcm(apply(denominator, Vec(f))), F = f*d^2);
  test3(F, Strprintf("RUNGBGP (%s, %s, %s)", s, t, u));
);
print("RUNGBGPDONE");
}
quit;
