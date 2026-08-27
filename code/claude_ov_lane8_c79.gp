/* Lane 8 (overnight 2026-07-25): Z/63 on the contact charts.
 *
 * A rational point of order 63 = 7*9 on J.  The two charts of this repo that
 * carry a marked rational class of order 7 resp. 9 on a quintic model with a
 * rational Weierstrass point at infinity are
 *
 *   contact-7 (2 params):  h7 = 1 - (7/2)x + a x^2 + b x^3 ,  f = (h7^2+(x-1)^7)/x^2
 *   contact-9 (1 param) :  h9 = 1 - (9/2)x + (63/8)x^2 - (105/16)x^3 + a x^4 ,
 *                          f  = (h9^2+(x-1)^9)/x^4
 *
 * (both monic quintics; the marked class is [P-infty] with P = (1, h(1)) ).
 * Since dim M_2 = 3, the contact-7 image is a divisor and the contact-9 image is
 * a curve, so their intersection -- the locus with BOTH a rational 7- and a
 * rational 9-class of Weierstrass-difference type -- has expected dimension 0.
 *
 * Stage 1 here is the local landscape: for each good prime p, how many chart
 * points over F_p satisfy the necessary condition 63 | #J(F_p)?  A prime with
 * ZERO would close the chart outright (the recorded precedent: [35] dies at p=3
 * on contact-7+5).
 */

c9poly(a) = {
  my(h, g, f);
  h = 1 - 9/2*x + 63/8*x^2 - 105/16*x^3 + a*x^4;
  g = h^2 + (x-1)^9;
  f = sum(i=4, 9, polcoeff(g,i)*x^(i-4));
  f;
}

c7poly(a, b) = {
  my(h, g, f);
  h = 1 - 7/2*x + a*x^2 + b*x^3;
  g = h^2 + (x-1)^7;
  f = sum(i=2, 7, polcoeff(g,i)*x^(i-2));
  f;
}

\\ #J(F_p) for a quintic f with t_INTMOD coefficients; 0 if singular/degenerate
nJ(f, p) = {
  my(L);
  if(poldegree(f) != 5, return(0));
  if(poldegree(gcd(f, deriv(f))) > 0, return(0));
  L = hyperellcharpoly(f);
  subst(L, x, 1);
}

\\ contact-9: residues a mod p with 63 | #J(F_p); also report the marked-order sanity
c9scan(p) = {
  my(good = List(), nsm = 0, n9 = 0);
  for(A = 0, p-1,
    my(f = c9poly(Mod(A,p)), n);
    n = nJ(f, p);
    if(n == 0, next);
    nsm++;
    if(n % 9 == 0, n9++);
    if(n % 63 == 0, listput(good, A));
  );
  [nsm, n9, Vec(good)];
}

c7scan(p) = {
  my(good = List(), nsm = 0, n7 = 0);
  for(A = 0, p-1, for(B = 0, p-1,
    my(f = c7poly(Mod(A,p), Mod(B,p)), n);
    n = nJ(f, p);
    if(n == 0, next);
    nsm++;
    if(n % 7 == 0, n7++);
    if(n % 63 == 0, listput(good, [A,B]));
  ));
  [nsm, n7, Vec(good)];
}

c9landscape(plist) = {
  print("=== contact-9 chart (1 parameter a) : 63 | #J(F_p) landscape ===");
  for(i = 1, #plist,
    my(p = plist[i], r = c9scan(p));
    printf("p=%4d  smooth=%4d  9|#J: %4d  63|#J: %4d   good a: %s\n",
           p, r[1], r[2], #r[3], if(#r[3] <= 12, Str(r[3]), Str(#r[3], " values")));
  );
}

c7landscape(plist) = {
  print("=== contact-7 chart (2 parameters a,b) : 63 | #J(F_p) landscape ===");
  for(i = 1, #plist,
    my(p = plist[i], r = c7scan(p));
    printf("p=%4d  smooth=%6d  7|#J: %6d  63|#J: %6d  (density %.4f)\n",
           p, r[1], r[2], #r[3], 1.0*#r[3]/r[1]);
  );
}

/* ---- rational sieve on the 1-parameter contact-9 chart ---------------------
 * sweep a = n/d, |n| <= H, 1 <= d <= H, gcd(n,d)=1, and require 63 | #J(F_p)
 * for every sieving prime p not dividing d.  Survivors are printed for an exact
 * Magma TorsionSubgroup test.
 */
c9sieve(H, plist) = {
  my(bm = vector(#plist), cnt = 0, surv = List(), t0 = getwalltime());
  for(i = 1, #plist, bm[i] = Set(c9scan(plist[i])[3]));
  for(i = 1, #plist, printf("sieve prime %d: %d admissible residues\n", plist[i], #bm[i]));
  for(d = 1, H,
    for(n = -H, H,
      if(gcd(n,d) != 1, next);
      cnt++;
      my(ok = 1);
      for(i = 1, #plist,
        my(p = plist[i]);
        if(d % p == 0, next);
        if(!setsearch(bm[i], lift(Mod(n,p)/Mod(d,p))), ok = 0; break);
      );
      if(ok, listput(surv, [n,d]); printf("SURVIVOR a = %d/%d\n", n, d));
    );
    if(d % 50 == 0, printf("PROGRESS d=%d checked=%d survivors=%d wall=%ds\n",
                           d, cnt, #surv, getwalltime()-t0));
  );
  printf("SEARCH_DONE c9sieve H=%d checked=%d survivors=%d\n", H, cnt, #surv);
  Vec(surv);
}
