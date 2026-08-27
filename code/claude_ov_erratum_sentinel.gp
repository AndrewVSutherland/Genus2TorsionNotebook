/* claude_ov_erratum_sentinel.gp -- genus2red ERRATUM propagation (2026-07-25).
 *
 * Refines Lane 7's erratum.  genus2red does NOT unconditionally drop the
 * 2-part: it returns the sentinel exponent -1 at 2 in its factorisation matrix
 * when it cannot determine the conductor exponent at 2 (and then omits the
 * 2-part from N), but when it CAN determine it, it includes it.  Control curve
 * 9450.b.2 is an example of the latter (genus2red returns the full 9450, type
 * "[I{2-0-0}] page 170" at p = 2).
 *
 * The operational test is therefore: does the returned matrix contain a row
 * [2, -1]?  This script applies that test to every genus2red conductor
 * recorded in this repository: the eleven contact-7 three-root [2,2,14]
 * curves, and the BLP C4-corrected [2,22] curve of
 * notes/claude_generic_222_2214_plan_2026_07_23.md (recorded there as
 * "conductor odd part 645^2").
 *
 * Usage:   gp -q code/claude_ov_erratum_sentinel.gp
 * Markers: SENT / SENTSUMMARY / ERRATUM_SENTINEL_DONE
 *
 * GP TRAP (cost 10 minutes here, recording it): a function definition
 * "G7(v) = ..." must NOT go inside a {...} block.  Inside a block newlines are
 * ignored, so the function body swallows the ENTIRE rest of the block and the
 * script runs silently with no output and exit status 0.  Define functions
 * outside the block.
 */
G7(v) = -(v^5 - v^3 - v^2/2)/(v+1)^2;
{
trips = [[-10,-10/7,-1/2], [-5,-15/8,-15/22], [-3,-3/4,-3/5],
         [-15/8,-15/19,-1/2], [-5/18,-10/49,4/17], [-4/9,-4/25,4/17],
         [-511/61,-511/625,-1/2], [-165/41,-33/16,-165/289],
         [-164/297,-1/2,164/361], [-17/50,-34/189,34/121], [-1/2,-13/49,13/50]];
nsent = 0;
for(i=1, #trips,
  s = trips[i][1]; t = trips[i][2]; u = trips[i][3];
  c4 = (G7(s)-G7(t))/(s^2-t^2); c0 = G7(s) - c4*s^2;
  b = c4 - 2; a = 9/2 - c0 - c4;
  h = 1 - (7/2)*x + a*x^2 + b*x^3;
  f = (h^2 + (x-1)^7)/x^2;
  for(j=1, 3, if(subst(f, x, 1-trips[i][j]^2) != 0, error("root check failed")));
  /* clear denominators the CORRECT way (denominator() on a t_POL returns 1) */
  d = lcm(apply(denominator, Vec(f)));
  F = d^2 * f;
  G = genus2red(F);
  M = G[2];
  sent = 0;
  for(k=1, matsize(M)[1], if(M[k,1]==2 && M[k,2]==-1, sent = 1));
  if(sent, nsent++);
  print("SENT contact7 idx=", i, " (s,t,u)=(", s, ",", t, ",", u, ")  genus2red_N=", G[1],
        "  v2(genus2red_N)=", valuation(G[1],2), "  sentinel[2,-1]=", sent);
);
/* the BLP C4-corrected [2,22] curve */
fB = x^6-18*x^5-4001*x^4-22524*x^3+859039*x^2-1926258*x-9043839;
GB = genus2red(fB);
MB = GB[2]; sentB = 0;
for(k=1, matsize(MB)[1], if(MB[k,1]==2 && MB[k,2]==-1, sentB = 1));
print("SENT blpC4 [2,22]  genus2red_N=", GB[1], " = ", factor(GB[1]),
      "  is 645^2? ", GB[1]==645^2, "  sentinel[2,-1]=", sentB);
print("SENTSUMMARY  contact7 sentinel present in ", nsent, " of ", #trips,
      " ; blpC4 sentinel ", sentB);
print("ERRATUM_SENTINEL_DONE");
}
quit;
