/* Lane 8 (2026-07-25, third session): upgrade the transverse-deformation test
 * from LINES to polynomial ARCS.
 *
 * code/claude_ov_lane8_transverse.gp killed 156,695 straight lines through the
 * eight HLP/Howe anchor/target pairs.  The obvious gap is that a transverse
 * component of the torsion locus need not be a line: it could be a conic, a
 * cubic, ... .
 *
 * The rigorous criterion is unchanged and applies verbatim to ANY family
 * parametrised polynomially over Q:
 *
 *     f_t = f0 + t*g1 + t^2*g2 + ... + t^d*gd ,   f_0 = the anchor.
 *
 * If the order-N section extended over Q(t) along this arc, it would specialise
 * at EVERY t0 in F_p of good reduction, so N | #J_{t0}(F_p) for all such t0.
 * One failing residue kills the arc.  Reparametrising t does not change the
 * image, so degree-d arcs cover every rational curve through the anchor that
 * admits a polynomial parametrisation of degree <= d with the anchor at t = 0.
 *
 * Transversality is imposed on the LINEAR term g1: that is the tangent
 * direction of the arc at the anchor, i.e. exactly the condition that the arc
 * leaves the (2,2)-split locus to first order.  Arcs whose linear term is
 * TANGENT to the split locus but which leave it at higher order are counted
 * separately (the TANGENT-THEN-LEAVING row) -- those are precisely the arcs a
 * line test cannot see at all.
 *
 * usage: ANCHOR=k DEG=2 NARC=20000 gp -q -s 512M -f code/claude_ov_lane8_arcs.gp
 */

/* NB: gp's getenv returns the INTEGER 0 (not "") for an unset variable, so the
   naive  if(getenv(s)=="", default, eval(getenv(s)))  silently yields 0. */
env(s, def) = {my(v = getenv(s)); if(type(v) != "t_STR" || v == "", def, eval(v));}

ANC  = env("ANCHOR", 1);
DEG  = env("DEG",    2);
NARC = env("NARC",20000);
PK   = env("PKILL", 200);

{anchors = [
 ["Z/5xZ/10", 45600*x^6 + 289161*x^4 + 35186670*x^2 - 705688215, 50, [5], 1],
 ["Z/7xZ/7 ", x^6 + 3025*x^4 + 3232987*x^2 + 869675859,          49, [7], 1],
 ["Z/63    ", 897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496,   63, [],  1],
 ["Z/45    ", 13981*x^6 + 29240200*x^4 + 49996210000*x^2 + 168300000000, 45, [], 1],
 ["Z/35    ", 640*x^6 + 5040*x^4 + 2480*x^2 + 9295,              35, [],  1],
 ["Z/2xZ/24", x^6 + 46*x^4 + 409*x^2 + 840,                      48, [2], 1],
 ["Z/3xZ/12", -20*x^6 + 249*x^4 - 936*x^2 + 1296,                36, [3], 1],
 ["f70 Z/70", 3168*x^5 + 697*x^4 - 23220*x^3 + 37620*x^2 - 23328*x + 5184, 70, [], 0],
 ["f70 Z/35", 3168*x^5 + 697*x^4 - 23220*x^3 + 37620*x^2 - 23328*x + 5184, 35, [], 0]
];}

chiof(f, p, dg) = {my(g = Mod(1,p)*f); if(poldegree(g) != dg, return(0)); if(poldegree(gcd(g, deriv(g))) > 0, return(0)); hyperellcharpoly(g);}
compat(chi, N, ranks) = {my(nJ); if(chi == 0, return(-1)); nJ = subst(chi, x, 1); if(nJ % N != 0, return(0)); for(i = 1, #ranks, my(n = ranks[i]); if(nJ % (n^2) != 0, return(0)); if(((Mod(1,n)*chi) % (Mod(1,n)*(x-1)^2)) != 0, return(0))); 1;}

/* kill one arc: return [p,t] of the first residue proving it is not a
   Q(t)-family, or 0 if none was found below PK */
killarc(f0, gs, N, ranks, dg, D0, lc) = {
  forprime(p = 3, PK,
    if(D0 % p == 0 || lc % p == 0, next);
    for(t = 1, p-1,
      my(fp = f0 + sum(j = 1, #gs, t^j*gs[j]), c = chiof(fp, p, dg));
      if(c == 0, next);
      if(compat(c, N, ranks) == 0, return([p,t]))));
  0;
}

run() = {
  my(A = anchors[ANC], tag = A[1], f0 = A[2], N = A[3], ranks = A[4], even = A[5]);
  my(D0 = poldisc(f0), lc = pollead(f0), dg = poldegree(f0), nrm = 0);
  printf("############ ARCS of degree %d through anchor %d : %s   (N = %d)\n", DEG, ANC, tag, N);
  if(even,
    my(a = polcoeff(f0,6), b = polcoeff(f0,4), c = polcoeff(f0,2), d = polcoeff(f0,0));
    my(u1 = [2*c, 4*b, 6*a], u2 = [6*d, 4*c, 2*b]);
    nrm = [u1[2]*u2[3]-u1[3]*u2[2], u1[3]*u2[1]-u1[1]*u2[3], u1[1]*u2[2]-u1[2]*u2[1]];
    printf("  normal functional to the bielliptic locus = %s\n", nrm));
  /* sanity: the anchor itself must PASS at every good p < PK */
  my(bad0 = 0);
  forprime(p = 3, PK, if(D0 % p == 0 || lc % p == 0, next);
    if(compat(chiof(f0, p, dg), N, ranks) != 1, bad0++));
  printf("  sanity: anchor fails compat at %d good p < %d (must be 0)\n", bad0, PK);
  /* sanity: killarc must NOT kill the constant arc (g = 0 is excluded, so use
     the zero arc explicitly -- it must return 0 = not killed) */
  printf("  sanity: the CONSTANT arc f_t = f0 is not killed : %d\n",
         killarc(f0, vector(DEG, j, 0), N, ranks, dg, D0, lc) == 0);

  my(u1 = 0, u2 = 0);
  if(even,
    my(a = polcoeff(f0,6), b = polcoeff(f0,4), c = polcoeff(f0,2), d = polcoeff(f0,0));
    u1 = 2*c*x + 4*b*x^3 + 6*a*x^5;      /* = f0' */
    u2 = 6*d*x + 4*c*x^3 + 2*b*x^5);     /* = 6x f0 - x^2 f0' */

  my(ntr = 0, ktr = 0, ntg = 0, ktg = 0, worst = [0,0], surv = List());
  /* (A) arcs TRANSVERSE at first order */
  for(i = 1, NARC,
    my(gs = vector(DEG, j, sum(k = 0, dg, (random(41)-20)*x^k)));
    if(gs[1] == 0, next);
    if(even && polcoeff(gs[1],1)*nrm[1] + polcoeff(gs[1],3)*nrm[2] + polcoeff(gs[1],5)*nrm[3] == 0, next);
    ntr++;
    my(w = killarc(f0, gs, N, ranks, dg, D0, lc));
    if(w != 0,
      ktr++;
      if(w[1] > worst[1] || (w[1] == worst[1] && w[2] > worst[2]), worst = w),
      if(#surv < 5, listput(surv, gs)));
    if(i % 5000 == 0, printf("    PROGRESS transverse %d/%d killed %d\n", i, NARC, ktr)));
  /* (B) arcs TANGENT to the split locus at first order, leaving it at order 2.
     g1 is drawn from the FULL tangent space to the bielliptic locus:
     even part arbitrary, odd part in span(u1,u2).  A line test cannot see
     these -- the line f0 + t*g1 stays split to first order. */
  if(even && DEG >= 2,
    for(i = 1, NARC,
      my(e = sum(k = 0, 3, (random(41)-20)*x^(2*k)),
         r1 = random(41)-20, r2 = random(41)-20);
      my(g1 = e + r1*u1 + r2*u2);
      if(g1 == 0, next);
      my(gs = vector(DEG, j, if(j == 1, g1, sum(k = 0, dg, (random(41)-20)*x^k))));
      /* require the arc to actually leave the split locus at some later order */
      my(leaves = 0);
      for(j = 2, DEG,
        if(polcoeff(gs[j],1)*nrm[1] + polcoeff(gs[j],3)*nrm[2] + polcoeff(gs[j],5)*nrm[3] != 0,
           leaves = 1; break));
      if(!leaves, next);
      ntg++;
      my(w = killarc(f0, gs, N, ranks, dg, D0, lc));
      if(w != 0,
        ktg++;
        if(w[1] > worst[1] || (w[1] == worst[1] && w[2] > worst[2]), worst = w),
        if(#surv < 5, listput(surv, gs)));
      if(i % 5000 == 0, printf("    PROGRESS tangent %d/%d killed %d\n", i, NARC, ktg))));
  printf("  TRANSVERSE-AT-FIRST-ORDER arcs : %d tested, %d KILLED, %d survivors\n", ntr, ktr, ntr - ktr);
  printf("  TANGENT-THEN-LEAVING      arcs : %d tested, %d KILLED, %d survivors\n", ntg, ktg, ntg - ktg);
  printf("  hardest (p,t) needed = %s\n", worst);
  if(#surv, printf("  SURVIVING ARCS (first %d): %s\n", #surv, Vec(surv)));
  printf("SEARCH_DONE arcs anchor %d deg %d\n", ANC, DEG);
}
run();
quit
