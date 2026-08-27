\\ claude_ov_lane7_conductors.gp -- Lane 7 (overnight 2026-07-25).
\\
\\ Conductors of all ELEVEN contact-7 three-root curves, computed with PARI's
\\ genus2red (Liu's algorithm), fully INDEPENDENTLY of Magma.  Magma's
\\ Conductor() for these models emits
\\     "WARNING: Using Ogg's formula when v_2(D)>=12, no correctness guarantee"
\\ and returns a 2-part that genus2red does not confirm, so genus2red is the
\\ authority here (it is also what the 2026-07-23 note used).
\\
\\ Curves are rebuilt from scratch in gp from the lane-brief G-recipe.
\\ Markers: COND / CONDDONE

G(v) = -(v^5 - v^3 - v^2/2)/(v+1)^2;

trips = [ [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5], \
          [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17], \
          [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289], \
          [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], \
          [-1/2, -13/49, 13/50] ];

tags = ["known", "known", "known-RM", "known", "known", "known", \
        "NEW", "NEW", "NEW", "NEW", "NEW"];

{
for(i = 1, #trips,
  s = trips[i][1]; t = trips[i][2]; u = trips[i][3];
  c4 = (G(s) - G(t))/(s^2 - t^2);
  c0 = G(s) - c4*s^2;
  b  = c4 - 2;
  a  = 9/2 - c0 - c4;
  h  = 1 - 7/2*x + a*x^2 + b*x^3;
  num = h^2 + (x-1)^7;
  f  = num \ x^2;
  if (f*x^2 != num, error("not divisible by x^2"));
  \\ assert the three claimed rational Weierstrass points
  for(j = 1, 3, if (subst(f, x, 1 - trips[i][j]^2) != 0, error("root fail")));
  d  = denominator(f);
  F  = f * d^2;
  gr = genus2red(F);
  N  = gr[1];
  printf("COND %s (%s, %s, %s) N = %s = %s\n", tags[i], s, t, u, N, factor(N));
);
print("CONDDONE");
}
quit;
