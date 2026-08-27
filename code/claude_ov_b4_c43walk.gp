/* claude_ov_b4_c43walk.gp -- Lane 4 (route B4): DIRECT walk of the (4,3)
   component of Flynn's quadratic-factor incidence variety.
   Independent, member-by-member confirmation of everything the symbolic
   computation (code/claude_ov_b4_c43param.m) claims over Q(w):

     (a) factor type of F_{t(w)} is [1,2,3] for every non-degenerate member;
     (b) 2-rank of J(Q) is 1 for every such member;
     (c) NO member has a Galois-stable (2,2)-kernel: for every rational
         quadratic block, the resolvent cubic of the complementary quartic has
         no rational root (b4_resroots empty).  Equivalently the cubic factor
         has no rational root;
     (d) the raise test done directly in w-space (independent of the
         E = 92.a1 quotient sieve in claude_ov_b4_c43tables.gp): does the
         cubic  Dw*rho^3 - Dw*rho^2 + (Dw-Nw)*rho + Nw  have a rational root?
         Nw = a^3-2a^2b+2ab^2-b^3, Dw = a^2 b, w = a/b.

   Usage:
     echo 'HFULL=60;HRAISE=1200;read("code/claude_ov_b4_c43walk.gp")' | gp -q
*/

B4DIR = if(type(B4DIR)=="t_STR", B4DIR, "/home/claude/torsion_jac/code");
read(strprintf("%s/claude_ov_b4_richelot.gp", B4DIR));

HFULL  = if(type(HFULL)=="t_INT",  HFULL,  60);    /* box for the full check   */
HRAISE = if(type(HRAISE)=="t_INT", HRAISE, 1200);  /* box for the raise search */

x='x; rho='rho;

Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;
Tw(ww) = (ww^3-2*ww^2+2*ww-1)/ww^2;
Uw(ww) = (2*ww^3-2*ww^2+2*ww-1)/ww^2;
Vw(ww) = (ww^2-ww+1)^2/ww^2;
tw(ww) = -Tw(ww)^2;

/* ------------------------------------------------------------------ (a)-(c) */
{
  my(nm=0, nbadtype=0, nbadrank=0, nkernel=0, nblocks=0, nsing=0, types=Map());
  for(b = 1, HFULL,
    for(a = -HFULL, HFULL,
      if(a == 0 || gcd(a,b) != 1, next);
      my(ww = a/b, tt, f, degs, tr, bl, kern, key);
      if(ww == 1, next);                       /* t = 0, F = x^2(x^2+x+1)^2   */
      tt = tw(ww);
      f = Fl(tt);
      if(poldisc(f) == 0, nsing++; next);
      /* the parametrized quadratic must really divide f */
      if(f % (x^2 + Uw(ww)*x + Vw(ww)) != 0, error("block does not divide at w=", ww));
      nm++;
      degs = b4_factortype(f);
      key = Str(degs);
      mapput(types, key, if(mapisdefined(types,key), mapget(types,key), 0) + 1);
      if(degs != [1,2,3]~ && Vec(degs) != [1,2,3], nbadtype++;
         print("  TYPE w=", ww, " degs=", degs));
      tr = b4_tworank(degs);
      if(tr != 1, nbadrank++; print("  2RANK w=", ww, " degs=", degs, " tworank=", tr));
      bl = b4_blocks(f);
      nblocks += #bl;
      kern = 0;
      foreach(bl, B, if(#b4_resroots(B[2]) > 0, kern++));
      if(kern > 0, nkernel++;
         print("  KERNEL w=", ww, " blocks=", #bl, " with-resolvent-root=", kern));
    );
  );
  print("FULLCHECK box |a|,|b| <= ", HFULL, "  members=", nm, " singular-skipped=", nsing);
  print("  factor types seen: ", Vec(types), " -> ", [mapget(types,k) | k <- Vec(types)]);
  print("  members with factor type != [1,2,3]      : ", nbadtype);
  print("  members with 2-rank != 1                 : ", nbadrank);
  print("  rational quadratic blocks total          : ", nblocks);
  print("  members with a GALOIS-STABLE (2,2)-KERNEL: ", nkernel);
}

/* --------------------------------------------------------------------- (d) */
{
  my(nn=0, nhit=0, t0=getwalltime());
  for(b = 1, HRAISE,
    for(a = -HRAISE, HRAISE,
      if(a == 0 || gcd(a,b) != 1, next);
      my(Nw, Dw, P, R);
      Nw = a^3 - 2*a^2*b + 2*a*b^2 - b^3;
      Dw = a^2*b;
      nn++;
      P = Dw*rho^3 - Dw*rho^2 + (Dw-Nw)*rho + Nw;
      if(poldegree(P) < 1, next);
      /* cheap modular prefilter: a rational root survives mod every good p */
      my(ok = 1);
      foreach([7,11,13,17,19,23,29], pp,
        if(ok && Dw % pp != 0 && #polrootsmod(P, pp) == 0, ok = 0));
      if(!ok, next);
      R = select(z -> poldegree(z) == 1, Vec(factor(P)[,1]~));
      if(#R > 0,
        nhit++;
        print("  RAISE-CANDIDATE w=", a, "/", b, "  rho-cubic linear factors: ", R));
    );
  );
  print("RAISESEARCH box |a|,|b| <= ", HRAISE, "  w-values tested=", nn,
        "  hits=", nhit, "  wall=", (getwalltime()-t0)/1000.0, "s");
}
print("C43WALK_DONE");
quit;
