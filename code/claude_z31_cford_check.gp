/* Task B3 helper: exact over-Q CF order of D_inf for CAND lines from
 * claude_z31_cf31_scan sweep mode.  Same recursion/guards as the C core
 * (Platonov-Petrunin; total quasi-period degree INCLUDING deg a_0 = 3).
 * Usage:  gp -q claude_z31_cford_check.gp   (edit the fs vector below), or
 *         \r claude_z31_cford_check.gp  inside gp and call
 *         cford(x^6 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0, 40)
 * A CAND with cford(f,40) == 31 and issquarefree(f) is a genuine
 * ord(D_inf) = 31 curve; then hand it to Magma for TorsionSubgroup +
 * simplicity certification (validate-and-record-a-hit skill).
 */
cford(f, cap) = {
  my(s, Pi, Qi, tot, ai, Pn, Qn, d);
  s = x^3;
  for (k = 1, 3,
    d = f - s^2;
    if (poldegree(d) <= 2, break);
    s = s + (polcoeff(d, 6-k)/(2*polcoeff(s,3)))*x^(3-k));
  Pi = 0; Qi = 1; tot = 0;
  for (i = 0, cap+10,
    if (Qi == 0, return(0));
    ai = (Pi + s) \ Qi;
    if (poldegree(ai) < 0, return(0));
    if (i >= 1 && poldegree(ai) < 1, return(0));
    tot += poldegree(ai);
    if (tot > cap, return(0));
    Pn = ai*Qi - Pi;
    if ((f - Pn^2) % Qi != 0, return(0));
    Qn = (f - Pn^2) \ Qi;
    Pi = Pn; Qi = Qn;
    if (i >= 1 && poldegree(Qi) <= 0 && Qi != 0, return(tot)));
  return(0);
}
/* self-check on the validated vectors, then the two e2e sweep finds */
{
fs = [[x^6 + 6*x^4 + 4*x^3 + 9*x^2 + 4*x + 4, 14],
      [(x^2 - x + 1)*(x^4 - x^3 + 9*x^2 + 8*x - 8), 18],
      [x^6 + 2*x^5 - 5*x^4 - 14*x^3 - 3*x^2 + 24*x + 28, 7],
      [x^6 + 3/2*x^4 + 1/2*x^3 + 9/16*x^2 + 1/8*x + 1/16, 14],
      [x^6 + 3/2*x^4 + 5*x^3 + 9/16*x^2 + 31/4*x + 33/4, 14]];
for (i = 1, #fs,
  o = cford(fs[i][1], 200);
  printf("cford = %d expect %d sqfree = %d %s\n", o, fs[i][2],
         issquarefree(fs[i][1]), if(o == fs[i][2], "PASS", "FAIL")));
}
quit
