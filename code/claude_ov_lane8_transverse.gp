/* Lane 8 (2026-07-25, third session): the TRANSVERSE-DEFORMATION experiment.
 *
 * Lane item: "put the split point in a chart, write the torsion condition around
 * it, identify the split divisor through it, and determine whether the torsion
 * locus has a component TRANSVERSE to the split locus."
 *
 * ---------------------------------------------------------------------------
 * WHY FIRST-ORDER DEFORMATION THEORY IS VACUOUS HERE
 *
 * In characteristic 0 the group scheme J[N] is FINITE ETALE over the base.  A
 * rational point of order N on the anchor therefore lifts UNIQUELY over every
 * Artinian thickening and over the completed local ring of moduli at the anchor.
 * So the "torsion locus" is not a subvariety of M_2 at all: it is the image of
 * the rational points of the finite etale cover
 *      W_N = { (C,P) : P in J(C), ord(P) = N }  --->  M_2 ,
 * every component of which is 3-dimensional and unramified over M_2.  Hence:
 *   (i)  there is NO local obstruction to deforming off the split locus, and
 *   (ii) the tangent/normal computation carries no information (this repeats,
 *        for these anchors, exactly what notes/m612_hlp_deformation.md found
 *        for [6,12]: 21x28 Jacobian of rank 21, kernel 7 = 3 moduli + 4 group).
 * The whole question is arithmetic: does W_N have rational points off the
 * divisor pi^{-1}(split locus)?
 *
 * ---------------------------------------------------------------------------
 * WHAT *CAN* BE DECIDED, AND HOW
 *
 * Sharp, rigorous, cheap criterion.  Let g be a direction and consider the line
 *      C_t : y^2 = f0(x) + t*g(x) ,  t in Q.
 * If the order-N section extended over Q(t) (i.e. if this line lay inside the
 * torsion locus as a FAMILY), it would specialise at every t0 in F_p of good
 * reduction, so N | #J_{t0}(F_p) for EVERY such t0.  Therefore:
 *
 *     ONE residue t0 in F_p with good reduction and N nmid #J_{t0}(F_p)
 *     PROVES that the line f0 + t*g carries no Q(t)-rational order-N family.
 *
 * Positive control that the criterion has power: on the contact-7 chart (which
 * does carry a marked rational order-7 class) the density of 7 | #J(F_p) is
 * exactly 1.000 -- see results/claude_ov_lane8_c7_landscape.log and stage (0).
 *
 * Stages:
 *   (0) positive control on the contact-7 / contact-9 charts (density must be 1)
 *   (1) anchor sanity: N | #J(F_p) at every good p < PMAX
 *   (2) densities: 3 transverse coordinate lines  f0 + t*x^{1,3,5}
 *                  3 in-locus  coordinate lines  f0 + t*x^{0,2,4}   (stay bielliptic)
 *                  baseline over random sextics at the same primes
 *   (3) EXHAUSTIVE direction kill: every primitive sign-normalised odd direction
 *       g = c1*x + c3*x^3 + c5*x^5 with |ci| <= BOX that is transverse to the
 *       bielliptic locus, plus NRDIR random directions in the full 7-dim
 *       coefficient space.  For each, find a residue that kills the line.
 *
 * The tangent space to the bielliptic (= (2,2)-split) locus at an EVEN sextic
 * f0 = a x^6 + b x^4 + c x^2 + d is  {even sextics} + span(u1,u2)  where
 *      u1 = f0'             = (2c, 4b, 6a)   in coords (x, x^3, x^5)
 *      u2 = 6x f0 - x^2 f0' = (6d, 4c, 2b)
 * (the sl2 + scaling action on binary sextics).  It has codimension 1, so a
 * direction g is TRANSVERSE iff  oddpart(g) . (u1 x u2) != 0.
 *
 * usage: ANCHOR=k PMAX=300 BOX=12 gp -q -f code/claude_ov_lane8_transverse.gp
 *        ANCHOR=0 runs the control and all anchors; ANCHOR=-1 runs the control only.
 */

/* NB: gp's getenv returns the INTEGER 0 (not "") for an unset variable, so the
   naive  if(getenv(s)=="", default, eval(getenv(s)))  silently yields 0.  The
   committed runs all passed every variable explicitly, so they are unaffected;
   this makes the defaults actually work. */
env(s, def) = {my(v = getenv(s)); if(type(v) != "t_STR" || v == "", def, eval(v));}

ANC  = env("ANCHOR",  0);
PMAX = env("PMAX",  300);
BOX  = env("BOX",    12);
NRD  = env("NRDIR",3000);

/* ---- anchors: [tag, f0, N, rank-2 primes, bielliptic-even-model?] ---- */
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

/* chi(T) of y^2 = f over F_p.  dg > 0 forces NO degree drop (= good reduction of
   the given plane model); dg = 0 only asks for a smooth genus-2 fibre. */
chiof(f, p, dg) = {my(g = Mod(1,p)*f); if(dg, if(poldegree(g) != dg, return(0)), if(poldegree(g) < 5, return(0))); if(poldegree(gcd(g, deriv(g))) > 0, return(0)); hyperellcharpoly(g);}

/* torsion-compatibility: N | #J and, for n in ranks, n-rank >= 2 */
compat(chi, N, ranks) = {my(nJ); if(chi == 0, return(-1)); nJ = subst(chi, x, 1); if(nJ % N != 0, return(0)); for(i = 1, #ranks, my(n = ranks[i]); if(nJ % (n^2) != 0, return(0)); if(((Mod(1,n)*chi) % (Mod(1,n)*(x-1)^2)) != 0, return(0))); 1;}

/* ------------------------- stage (0) positive control ------------------- */
c7poly(a, b) = {my(h = 1 - 7/2*x + a*x^2 + b*x^3, g = h^2 + (x-1)^7); sum(i = 2, 7, polcoeff(g,i)*x^(i-2));}
c9poly(a) = {my(h = 1 - 9/2*x + 63/8*x^2 - 105/16*x^3 + a*x^4, g = h^2 + (x-1)^9); sum(i = 4, 9, polcoeff(g,i)*x^(i-4));}

control() = {
  printf("########## STAGE 0 : positive control -- charts that DO carry a marked class\n");
  for(pi = 1, 6,
    my(p = [11,13,17,19,23,29][pi], n7 = 0, t7 = 0, n9 = 0, t9 = 0);
    for(A = 0, p-1, for(B2 = 0, p-1,
      my(c = chiof(c7poly(Mod(A,p), Mod(B2,p)), p, 0));
      if(c == 0, next); t7++; if(subst(c,x,1) % 7 == 0, n7++)));
    for(A = 0, p-1,
      my(c = chiof(c9poly(Mod(A,p)), p, 0));
      if(c == 0, next); t9++; if(subst(c,x,1) % 9 == 0, n9++));
    printf("  p=%3d  contact-7: 7|#J at %d/%d = %.4f   contact-9: 9|#J at %d/%d = %.4f\n",
           p, n7, t7, if(t7, n7/t7*1., 0), n9, t9, if(t9, n9/t9*1., 0)));
  printf("  => a line inside the torsion locus gives density EXACTLY 1; the test has power.\n\n");
}

/* --------------------------------- main -------------------------------- */
run(k) = {
  my(A = anchors[k], tag = A[1], f0 = A[2], N = A[3], ranks = A[4], even = A[5]);
  my(D0 = poldisc(f0), lc = pollead(f0), dg = poldegree(f0), nrm = 0, u1, u2);
  printf("############ anchor %d : %s\n  f0 = %s\n", k, tag, f0);
  printf("  N = %d   rank-2 primes = %s   even/bielliptic model = %d   deg f0 = %d\n", N, ranks, even, dg);
  if(even,
    my(a = polcoeff(f0,6), b = polcoeff(f0,4), c = polcoeff(f0,2), d = polcoeff(f0,0));
    u1 = [2*c, 4*b, 6*a]; u2 = [6*d, 4*c, 2*b];
    nrm = [u1[2]*u2[3]-u1[3]*u2[2], u1[3]*u2[1]-u1[1]*u2[3], u1[1]*u2[2]-u1[2]*u2[1]];
    printf("  bielliptic tangent (odd part): u1 = %s  u2 = %s\n  NORMAL functional = %s\n", u1, u2, nrm));

  /* ---- stage 1: anchor sanity ---- */
  my(ok = 1, ngood = 0);
  forprime(p = 3, PMAX, if(D0 % p == 0 || lc % p == 0, next);
    my(c = chiof(f0, p, dg)); ngood++;
    if(compat(c, N, ranks) != 1, ok = 0; printf("  !! anchor FAILS compat at p=%d\n", p)));
  printf("  STAGE1 anchor passes the compatibility test at all %d good p < %d : %d\n", ngood, PMAX, ok);

  /* ---- stage 2: densities ---- */
  my(hitL = 0, totL = 0, hitF = 0, totF = 0, hitR = 0, totR = 0, firstfail = vector(6));
  forprime(p = 3, PMAX, if(D0 % p == 0 || lc % p == 0, next);
    for(j = 1, 6,
      my(kk = [1,3,5,0,2,4][j]);
      for(t = 1, p-1,
        my(c = chiof(f0 + t*x^kk, p, dg), r);
        if(c == 0, next); r = compat(c, N, ranks);
        if(j <= 3, totL++; if(r, hitL++), totF++; if(r, hitF++));
        if(!r && firstfail[j] == 0, firstfail[j] = [p, t])));
    for(i = 1, 4*p,
      my(f = sum(j = 0, 6, random(p)*x^j), c = chiof(f, p, 0));
      if(c == 0, next); totR++; if(compat(c, N, ranks), hitR++));
    if(p % 47 < 8, printf("    PROGRESS p=%d  lines %d/%d  flat %d/%d  rand %d/%d\n",
                          p, hitL, totL, hitF, totF, hitR, totR)));
  printf("  STAGE2 TRANSVERSE lines f0+t*x^{1,3,5} : %d/%d = %.5f\n", hitL, totL, if(totL, hitL/totL*1., 0));
  printf("  STAGE2 IN-LOCUS   lines f0+t*x^{0,2,4} : %d/%d = %.5f\n", hitF, totF, if(totF, hitF/totF*1., 0));
  printf("  STAGE2 BASELINE   random sextics       : %d/%d = %.5f\n", hitR, totR, if(totR, hitR/totR*1., 0));
  for(j = 1, 6, printf("  STAGE2 first killing residue (p,t) for f0+t*x^%d : %s\n", [1,3,5,0,2,4][j], firstfail[j]));

  /* ---- stage 3: exhaustive odd-direction kill (bielliptic anchors only) ---- */
  if(even,
    my(ndir = 0, nkill = 0, worst = [0,0], surv = List());
    for(c1 = -BOX, BOX, for(c3 = -BOX, BOX, for(c5 = -BOX, BOX,
      if(c1 == 0 && c3 == 0 && c5 == 0, next);
      if(gcd(gcd(c1,c3),c5) != 1, next);
      if(c1 < 0 || (c1 == 0 && (c3 < 0 || (c3 == 0 && c5 < 0))), next);
      if(c1*nrm[1] + c3*nrm[2] + c5*nrm[3] == 0, next);
      ndir++;
      my(gg = c1*x + c3*x^3 + c5*x^5, done = 0);
      forprime(p = 3, 200, if(D0 % p == 0 || lc % p == 0, next);
        for(t = 1, p-1,
          my(c = chiof(f0 + t*gg, p, dg));
          if(c == 0, next);
          if(!compat(c, N, ranks), done = 1;
             if(p > worst[1] || (p == worst[1] && t > worst[2]), worst = [p,t]); break));
        if(done, break));
      if(done, nkill++, listput(surv, gg)))));
    printf("  STAGE3 exhaustive odd box |ci|<=%d : %d transverse primitive directions, %d KILLED, %d survivors\n",
           BOX, ndir, nkill, #surv);
    printf("  STAGE3 hardest (p,t) needed = %s\n", worst);
    if(#surv, printf("  STAGE3 SURVIVING DIRECTIONS: %s\n", Vec(surv))));

  /* ---- stage 3R: random directions in the full coefficient space ---- */
  my(nr = 0, nrk = 0, rsurv = List(), rworst = [0,0]);
  for(i = 1, NRD,
    my(gg = sum(j = 0, dg, (random(41)-20)*x^j));
    if(gg == 0, next);
    if(even && (polcoeff(gg,1)*nrm[1] + polcoeff(gg,3)*nrm[2] + polcoeff(gg,5)*nrm[3]) == 0, next);
    nr++;
    my(done = 0);
    forprime(p = 3, 200, if(D0 % p == 0 || lc % p == 0, next);
      for(t = 1, p-1,
        my(c = chiof(f0 + t*gg, p, dg));
        if(c == 0, next);
        if(!compat(c, N, ranks), done = 1;
           if(p > rworst[1] || (p == rworst[1] && t > rworst[2]), rworst = [p,t]); break));
      if(done, break));
    if(done, nrk++, listput(rsurv, gg)));
  printf("  STAGE3R random directions, coeffs in [-20,20]^%d : %d transverse, %d KILLED, %d survivors\n",
         dg+1, nr, nrk, #rsurv);
  printf("  STAGE3R hardest (p,t) needed = %s\n", rworst);
  if(#rsurv, printf("  STAGE3R SURVIVING DIRECTIONS: %s\n", Vec(rsurv)));
  printf("SEARCH_DONE anchor %d %s\n", k, tag);
}

if(ANC == 0, control(); for(k = 1, #anchors, run(k)), if(ANC == -1, control(), run(ANC)));
quit
