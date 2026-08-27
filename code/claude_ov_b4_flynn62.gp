/* claude_ov_b4_flynn62.gp -- Lane 4 (route B4) seed stream:
   enumerate members of the (6,2) component of the quadratic-factor incidence
   of Flynn's order-11 family, and test the Galois-stable-Richelot condition
   "resolvent cubic of the quartic cofactor has a rational root".

   Flynn family (J. Number Theory 36 (1990) Application 3.1), repo form:
     F_t(x) = x^6+2x^5+(2t+3)x^4+2x^3+(t^2+1)x^2+2t(1-t)x+t^2
   Incidence {x^2+u x+v | F_t}, (6,2) component:
     a2(u) v^2 + b1(u) v + c0(u) = 0,
     a2 = (u-2)^4, b1 = 2(u^5-5u^4+9u^3-8u^2+4u-2), c0 = u^2(u^2-u+1)^2.
   Its smooth model is E: Y^2 = w^3+3w^2+2w+1 (cond 92, rank 1) with u = -w,
   sqrt(b1^2-4 a2 c0) = 4(u-1)Y.

   Usage:  gp -q -D parisize=4000000000 claude_ov_b4_flynn62.gp
   Optional: set NMAX below (default 30).
*/

NMAX = if(type(NMAX)=="t_INT", NMAX, 30);
OUTF = "/home/claude/torsion_jac/results/claude_ov_b4_flynn62_rows.txt";

x='x; t='t; u='u; v='v; z='z;

Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;
a2(uu) = (uu-2)^4;
b1(uu) = 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);
c0(uu) = uu^2*(uu^2-uu+1)^2;

/* ---- 0. verify the discriminant identity symbolically (PARI factor() drops
        the content, so re-attach it explicitly and compare polynomials) ---- */
{
  D = b1(u)^2 - 4*a2(u)*c0(u);
  Dtarget = 16*(u-1)^2*subst(w^3+3*w^2+2*w+1, w, -u);
  print("DISCCHECK  b1^2-4a2c0 - 16(u-1)^2*E(-u) = ", D - Dtarget);
}

/* ---- 1. resolvent cubic of a monic quartic; rational-root test ---- */
resolvent(q4) = {
  my(a3,a2c,a1,a0,p,q,r);
  a3 = polcoef(q4,3,x); a2c = polcoef(q4,2,x);
  a1 = polcoef(q4,1,x); a0 = polcoef(q4,0,x);
  p = a2c - 3*a3^2/8;
  q = a1 - a3*a2c/2 + a3^3/8;
  r = a0 - a1*a3/4 + a2c*a3^2/16 - 3*a3^4/256;
  z^3 - p*z^2 - 4*r*z + (4*p*r - q^2);
};

hasratroot(c) = {
  my(f);
  if(poldegree(c) < 1, return(0));
  f = factor(c);
  for(i=1, matsize(f)[1], if(poldegree(f[i,1]) == 1, return(1)));
  0;
};

/* ---- 2. positive control for the whole D4/Richelot logic ---- */
{
  my(qc);
  qc = x^4-4*x^2+2;                       /* the verified Magma control */
  print("CONTROL x^4-4x^2+2 : Gal<=D4 by resolvent? ", hasratroot(resolvent(qc)),
        "   polgalois=", polgalois(qc));
  qc = x^4+x^3+x^2+x+1;                   /* C4 (cyclotomic), should pass */
  print("CONTROL x^4+x^3+x^2+x+1 : ", hasratroot(resolvent(qc)),
        "   polgalois=", polgalois(qc));
  qc = x^4+x+1;                           /* S4, should fail */
  print("CONTROL x^4+x+1 : ", hasratroot(resolvent(qc)),
        "   polgalois=", polgalois(qc));
}

/* ---- 3. RM prefilter: squarefree cores of the real Weil subfield disc ---- */
discsig(f6, np) = {
  my(S=List(), p=3, n=0, cp, c3, c2, d);
  while(n < np && p < 400,
    if(poldisc(f6) != 0 && Mod(poldisc(f6),p) != 0 && poldegree(f6)>=5,
      cp = hyperellcharpoly(Mod(f6,p));
      c3 = polcoef(cp,3); c2 = polcoef(cp,2);
      d = c3^2 - 4*(c2 - 2*p);
      if(d != 0, listput(S, core(d)); n++);
    );
    p = nextprime(p+1);
  );
  Set(Vec(S));
};

/* ---- 4. enumerate the stream ---- */
{
  my(E, G, seen, nmem, npass, rows);
  E = ellinit([0,3,0,2,1]);
  print("E = ", E.disc, "  conductor=", ellglobalred(E)[1],
        "  torsion=", elltors(E)[1]);
  G = [0,1];
  print("G=(0,1) on E? ", ellisoncurve(E,G), "  order=", ellorder(E,G),
        "  canonical height=", ellheight(E,G));

  seen = Map(); nmem = 0; npass = 0; rows = List();
  for(n = -NMAX, NMAX,
    if(n == 0, next);
    my(P, w, Y, uu);
    P = ellmul(E, G, n);
    if(#P < 2, next);
    w = P[1]; Y = P[2];
    uu = -w;
    if(uu == 2 || uu == 1, next);
    for(sg = 0, 1,
      my(A,B,vv,ff,dv,r,ee1,ee0,g,rts);
      A = a2(uu); B = b1(uu);
      vv = (-B + (-1)^sg*4*(uu-1)*Y)/(2*A);
      if(vv == 0, next);
      ff = Fl(t);
      dv = divrem(ff, x^2+uu*x+vv);
      r = dv[2];
      ee1 = polcoef(r,1,x); ee0 = polcoef(r,0,x);
      if(ee1 == 0 && ee0 == 0, next);
      g = if(ee1==0, ee0, if(ee0==0, ee1, gcd(ee1,ee0)));
      if(type(g) != "t_POL" || poldegree(g,t) < 1, next);
      rts = nfroots(, g);
      foreach(rts, t0,
        my(key, F, dq, q4, fa, degs, res, ok, sig);
        key = Str([uu,vv,t0]);
        if(mapisdefined(seen,key), next);
        mapput(seen,key,1);
        F = subst(ff, t, t0);
        if(poldegree(F) < 5 || poldisc(F) == 0, next);
        dq = divrem(F, x^2+uu*x+vv);
        if(dq[2] != 0, print("NODIV n=",n," u=",uu); next);
        q4 = dq[1];
        fa = factor(q4);
        degs = vecsort(apply(pp->poldegree(pp), fa[,1]~));
        nmem++;
        res = resolvent(q4);
        ok = if(degs == [4], hasratroot(res), -1);
        if(ok == 1, npass++);
        print("MEMBER n=", n, " sg=", sg,
              " hgt=", #Str(numerator(uu))+#Str(denominator(uu)),
              " q4type=", degs,
              " D4=", ok,
              " u=", uu);
        listput(rows, [n, sg, uu, vv, t0, degs, ok]);
      );
    );
  );
  print("TOTALMEMBERS ", nmem, "  D4PASS ", npass);
  write(OUTF, "");
  foreach(Vec(rows), rw, write(OUTF, rw));
  print("SEARCH_DONE flynn62 NMAX=", NMAX);
}
quit
