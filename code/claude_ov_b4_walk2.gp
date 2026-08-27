/* claude_ov_b4_walk2.gp -- Lane 4 (route B4): DIRECT Richelot walk of the
   Flynn (6,2) stream, sharded.  For every member it
     - builds F_t and the rational quadratic block,
     - applies the Galois-stable Richelot isogeny over Q,
     - FACTORS the codomain sextic and records its 2-rank,
     - and independently evaluates the symbolic predictors A1, D1
       (code/claude_ov_b4_locus_data.gp) so the two agree or we learn something.

   Members <-> points of E(Q)\{O}: for n >= 1 and sg in {0,1},
       u = -x(nG),  yy = (-1)^sg*y(nG),  v = (-b1(u) + 4(u-1)yy)/(2 a2(u)).
   Sweeping n = 1..NHI with both sg enumerates every member with |n| <= NHI.

   Usage: echo 'NLO=1;NHI=150;NOFF=0;NSTEP=12;TAG="w0";read("code/claude_ov_b4_walk2.gp")' | gp -q
*/

B4DIR = if(type(B4DIR)=="t_STR", B4DIR, "/home/claude/torsion_jac/code");
read(strprintf("%s/claude_ov_b4_richelot.gp", B4DIR));
read(strprintf("%s/claude_ov_b4_locus_data.gp", B4DIR));

NLO  = if(type(NLO)=="t_INT", NLO, 1);
NHI  = if(type(NHI)=="t_INT", NHI, 80);
NOFF = if(type(NOFF)=="t_INT", NOFF, 0);
NSTEP= if(type(NSTEP)=="t_INT", NSTEP, 1);

x='x; t='t; u='u;
Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;
a2c(uu) = (uu-2)^4;
b1c(uu) = 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);

dlcm(P) = { my(L=1); for(i=0, poldegree(P), L = lcm(L, denominator(polcoef(P,i,u)))); L; };
mkloc(C0N, C1N, DEN) = { my(L,M); L = lcm(dlcm(C0N), dlcm(C1N)); M = dlcm(DEN);
                         [L*C0N, L*C1N, M*DEN, L*M]; };
LOCA = mkloc(A1_C0_NUM, A1_C1_NUM, A1_C0_DEN);
LOCD = mkloc(D1_C0_NUM, D1_C1_NUM, D1_C0_DEN);
locval(LOC, uu, yyv) = LOC[4]*(subst(LOC[1],u,uu) + subst(LOC[2],u,uu)*yyv)*subst(LOC[3],u,uu);

{
  my(E, G, nm=0, nraise=0, nfail=0, nmis=0);
  E = ellinit([0,3,0,2,1]); G = [0,1];
  for(n = NLO, NHI,
    if(n % NSTEP != NOFF, next);
    if(n == 0, next);
    my(P, uu, Y0);
    P = ellmul(E, G, n);
    if(#P < 2, next);
    uu = -P[1]; Y0 = P[2];
    if(uu == 2 || uu == 1, next);
    for(sg = 0, 1,
      my(yy, A, B, vv, ff, dv, r, ee1, ee0, g, rts, F, bl, rr, out, degs, tr,
         a1v, d1v, sqA1, sqD1, pred, obs, hgt, deg2);
      yy = (-1)^sg*Y0;
      A = a2c(uu); B = b1c(uu);
      vv = (-B + 4*(uu-1)*yy)/(2*A);
      if(vv == 0, next);
      ff = Fl(t);
      dv = divrem(ff, x^2+uu*x+vv);
      r = dv[2];
      ee1 = polcoef(r,1,x); ee0 = polcoef(r,0,x);
      if(ee1 == 0 && ee0 == 0, next);
      g = if(ee1==0, ee0, if(ee0==0, ee1, gcd(ee1,ee0)));
      if(type(g) != "t_POL" || poldegree(g,t) < 1, next);
      rts = nfroots(, g);
      if(#rts == 0, next);
      F = subst(ff, t, rts[1]);
      if(poldegree(F) < 5 || poldisc(F) == 0, next);
      nm++;
      hgt = #Str(numerator(uu)) + #Str(denominator(uu));
      a1v = locval(LOCA, uu, yy); d1v = locval(LOCD, uu, yy);
      sqA1 = if(a1v == 0, -1, issquare(a1v));
      sqD1 = if(d1v == 0, -1, issquare(d1v));
      pred = if(sqA1==1, "RAISE", if(sqD1==1, "[1,1,4]", if(sqA1<0||sqD1<0, "DEGEN", "[2,4]")));
      bl = b4_blocks(F);
      if(#bl == 0, print("NOBLOCK n=",n," sg=",sg); nfail++; next);
      rr = b4_resroots(bl[1][2]);
      if(#rr == 0, print("NOD4 n=",n," sg=",sg); nfail++; next);
      out = b4_richelot(F, bl[1][1], bl[1][2], 1);
      if(out[1] != 1, print("DEGEN n=",n," sg=",sg," code=",out[1]); nfail++; next);
      degs = b4_factortype(out[2]);
      tr = b4_tworank(degs);
      obs = if(tr >= 2, "RAISE", Str(degs));
      if(tr >= 2,
        nraise++;
        print("RAISE n=", n, " sg=", sg, " codtype=", degs, " 2rank=", tr);
        print("RAISE_F ", b4_intmodel(out[2])));
      if(pred != obs && !(pred=="[2,4]" && obs=="[2, 4]") && !(pred=="[1,1,4]" && obs=="[1, 1, 4]") && pred != "DEGEN",
        nmis++; print("MISMATCH n=",n," sg=",sg," pred=",pred," obs=",obs));
      print("MEM n=", n, " sg=", sg, " uhgt=", hgt, " codtype=", degs,
            " cod2rank=", tr, " pred=", pred);
    );
    if(n % 10 == 0, print("PROGRESS n=", n, " members=", nm, " raises=", nraise));
  );
  print("WALK_DONE off=", NOFF, "/", NSTEP, " range=[", NLO, ",", NHI, "] members=", nm,
        " raises=", nraise, " mismatches=", nmis, " failures=", nfail);
}
quit
