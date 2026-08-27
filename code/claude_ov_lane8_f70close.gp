/* Lane 8 (2026-07-25, third session): CLOSING the f70 thread.
 *
 * db8137a established, at 299 good primes p < 2000, that
 *      #J(F_p) = #E1(F_p) * #E2(F_p)
 * for Howe's order-70 curve f70 and E1 = 66.c ([1,0,0,-45,81], torsion Z/10),
 * E2 = 858.k ([1,0,0,-5774401,5346023177], torsion Z/7); and that the gluing is
 * (3,3) rather than (2,2).  Here we upgrade #J(F_p) (one value) to the FULL
 * L-polynomial identity
 *      chi_J(T) = (T^2 - a_p(E1) T + p) (T^2 - a_p(E2) T + p)
 * at every good p below PMAX, which is the complete local statement of
 * J ~ E1 x E2 at p, and we re-run the mod-l congruence scan on a big range.
 *
 * usage: PMAX=20000 gp -q -f code/claude_ov_lane8_f70close.gp
 */

PMAX = if(getenv("PMAX")=="", 20000, eval(getenv("PMAX")));

f  = 3168*x^5 + 697*x^4 - 23220*x^3 + 37620*x^2 - 23328*x + 5184;   /* 144*f70 */
E1 = ellinit([1,0,0,-45,81]);
E2 = ellinit([1,0,0,-5774401,5346023177]);
bad = [2,3,11,13];

printf("f70 integral model : %s\n", f);
printf("E1 = [1,0,0,-45,81]                   conductor %d  torsion %s\n", ellglobalred(E1)[1], elltors(E1)[2]);
printf("E2 = [1,0,0,-5774401,5346023177]      conductor %d  torsion %s\n", ellglobalred(E2)[1], elltors(E2)[2]);
printf("bad primes: %s   PMAX = %d\n\n", bad, PMAX);

{
my(np = 0, agree = 0, disagree = List(), lmax = 40, cong = vector(lmax,i,1), congn = vector(lmax,i,0));
forprime(p = 5, PMAX,
  if(setsearch(Set(bad), p), next);
  my(g = Mod(1,p)*f);
  if(poldegree(g) != 5 || poldegree(gcd(g, deriv(g))) > 0, next);
  my(chi = hyperellcharpoly(g), a1 = ellap(E1,p), a2 = ellap(E2,p));
  my(prd = (x^2 - a1*x + p)*(x^2 - a2*x + p));
  np++;
  if(chi == prd, agree++, if(#disagree < 8, listput(disagree, [p, chi, prd])));
  for(l = 2, lmax, if(isprime(l), congn[l]++; if((a1 - a2) % l != 0, cong[l] = 0)));
  if(p % 5000 < 20, printf("  PROGRESS p=%d  tested %d  agree %d\n", p, np, agree));
);
printf("\nFULL L-POLYNOMIAL IDENTITY chi_J = L(E1)*L(E2) at good p < %d :\n", PMAX);
printf("  primes tested = %d   AGREE = %d   DISAGREE = %d\n", np, agree, np - agree);
if(#disagree, printf("  first disagreements: %s\n", Vec(disagree)));
printf("\nmod-l congruence a_p(E1) = a_p(E2) mod l for ALL good p < %d :\n", PMAX);
for(l = 2, lmax, if(isprime(l), printf("  l=%2d : %s  (%d primes tested)\n", l, if(cong[l], "HOLDS", "fails"), congn[l])));
}
printf("SEARCH_DONE f70close\n");
quit
