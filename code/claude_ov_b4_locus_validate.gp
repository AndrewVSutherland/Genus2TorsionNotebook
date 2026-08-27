/* claude_ov_b4_locus_validate.gp -- Lane 4 (route B4).
   VALIDATE the symbolically-derived raise-locus functions A1, D1 (dumped by
   code/claude_ov_b4_dump.m, extracted into code/claude_ov_b4_locus_data.gp)
   against the DIRECT Richelot walk on the (6,2) stream.

   Member <-> point of E(Q), E: Y^2 = w^3+3w^2+2w+1 (92.a1, rank 1, tors 1),
   E(Q) = Z.G, G=(0,1).  For n in Z\{0} and sg in {0,1}:
       u = -x(nG),  yy = (-1)^sg * y(nG),   yy^2 = -u^3+3u^2-2u+1.
   (sg=1 at n is the same member as sg=0 at -n; the walk enumerated both, so
   the honest indexing is by the point +-nG, i.e. by n in Z\{0}.)

   PREDICTION (from the symbolic derivation over K=Q(E)):
     codomain sextic = (rational quadratic h1) * (quartic q4'), generically
     factor type [2,4], 2-rank 1;
     q4' splits into two rational quadratics  <=>  A1 is a rational square;
     h1 splits into two rational linears      <=>  D1 is a rational square.
   OBSERVED: run the actual Richelot map and factor the codomain.

   Usage:  echo 'NLO=-40;NHI=-1;read("code/claude_ov_b4_locus_validate.gp")' | gp -q
*/

B4DIR = if(type(B4DIR)=="t_STR", B4DIR, "/home/claude/torsion_jac/code");
read(strprintf("%s/claude_ov_b4_richelot.gp", B4DIR));
read(strprintf("%s/claude_ov_b4_locus_data.gp", B4DIR));

NLO = if(type(NLO)=="t_INT", NLO, -40);
NHI = if(type(NHI)=="t_INT", NHI, -1);
DOWALK = if(type(DOWALK)=="t_INT", DOWALK, 1);

x='x; t='t;
Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;
a2c(uu) = (uu-2)^4;
b1c(uu) = 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);

evrat(nu, de, uu) = { my(d = subst(de,u,uu)); if(d==0, return(0)); subst(nu,u,uu)/d; };

A1val(uu, yy) = evrat(A1_C0_NUM,A1_C0_DEN,uu) + evrat(A1_C1_NUM,A1_C1_DEN,uu)*yy;
D1val(uu, yy) = evrat(D1_C0_NUM,D1_C0_DEN,uu) + evrat(D1_C1_NUM,D1_C1_DEN,uu)*yy;
Aval (uu, yy) = evrat(A_C0_NUM ,A_C0_DEN ,uu) + evrat(A_C1_NUM ,A_C1_DEN ,uu)*yy;

{
  my(E, G, nm=0, nmis=0, nraise=0);
  E = ellinit([0,3,0,2,1]);
  G = [0,1];
  print("E disc=", E.disc, " cond=", ellglobalred(E)[1], " tors=", elltors(E)[1]);
  for(n = NLO, NHI,
    if(n == 0, next);
    my(P, w0, Y0, uu);
    P = ellmul(E, G, n);
    if(#P < 2, next);
    w0 = P[1]; Y0 = P[2];
    uu = -w0;
    if(uu == 2 || uu == 1, next);
    for(sg = 0, 1,
      my(yy, A, B, vv, ff, dv, r, ee1, ee0, g, rts, a1v, d1v, sqA1, sqD1, pred);
      yy = (-1)^sg * Y0;
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
      nm++;
      a1v = A1val(uu,yy); d1v = D1val(uu,yy);
      sqA1 = issquare(a1v); sqD1 = issquare(d1v);
      pred = if(sqA1, "RAISE", if(sqD1, "[1,1,4]", "[2,4]"));
      if(sqA1, nraise++);
      if(DOWALK,
        my(F, bl, rr, out, degs, tr, obs, hgt);
        F = subst(ff, t, rts[1]);
        bl = b4_blocks(F);
        if(#bl == 0, print("NOBLOCK n=",n," sg=",sg); next);
        rr = b4_resroots(bl[1][2]);
        if(#rr == 0, print("NOD4 n=",n," sg=",sg); next);
        out = b4_richelot(F, bl[1][1], bl[1][2], 1);
        if(out[1] != 1, print("DEGEN n=",n," sg=",sg," code=",out[1]); next);
        degs = b4_factortype(out[2]);
        tr = b4_tworank(degs);
        obs = if(tr >= 2, "RAISE", if(degs == [1,1,4], "[1,1,4]", if(degs == [2,4], "[2,4]", Str(degs))));
        hgt = #Str(numerator(uu)) + #Str(denominator(uu));
        if(obs != pred,
          nmis++;
          print("MISMATCH n=", n, " sg=", sg, " pred=", pred, " obs=", obs,
                " degs=", degs, " sqA1=", sqA1, " sqD1=", sqD1),
          print("OK n=", n, " sg=", sg, " type=", obs, " uhgt=", hgt,
                " sqA1=", sqA1, " sqD1=", sqD1));
      ,
        print("PRED n=", n, " sg=", sg, " pred=", pred, " sqA1=", sqA1, " sqD1=", sqD1);
      );
    );
  );
  print("VALIDATE_DONE members=", nm, " mismatches=", nmis, " predicted_raises=", nraise,
        " range=[", NLO, ",", NHI, "]");
}
quit
