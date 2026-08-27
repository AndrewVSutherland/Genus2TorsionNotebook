/* claude_ov_b4_mktables.gp -- Lane 4 (route B4): build the mod-p Legendre
   tables for the three raise-locus double covers of E = 92.a1.

   E : Y^2 = w^3+3w^2+2w+1,  E(Q) = Z.G, G = (0,1), rank 1, trivial torsion.
   Member of the (6,2) stream at n in Z\{0}:  u = -x(nG), yy = y(nG).
   ((u,yy) and (u,-yy) are the sg=0/sg=1 pair, i.e. nG and -nG, so sweeping
    n over Z\{0} enumerates EVERY member exactly once.)

   Raise-locus test functions (square <=> the corresponding cover has a point):
     SQA1 = cA * (p0A(u) + p1A(u)*yy) * dA(u)      [q4' splits into 2 quadratics]
     SQD1 = cD * (p0D(u) + p1D(u)*yy) * dD(u)      [codomain quadratic splits]
     SQW  = cW * (p0W(u) + p1W(u)*yy) * dW(u)      [disc(q4') square]   (optional)

   Output (results/claude_ov_b4_tab_<tag>.txt), one block per prime:
     p N_p
     <N_p characters, '1' = Legendre symbol +1 or 0 (ALLOWED), '0' = -1 (KILLED)>
   with the r-th character corresponding to n = r (mod N_p).

   Usage:
     echo 'PMAX=600;TAG="small";read("code/claude_ov_b4_mktables.gp")' | gp -q
     echo 'PLIST=[...];TAG="mw";read("code/claude_ov_b4_mktables.gp")'  | gp -q
*/

B4DIR = if(type(B4DIR)=="t_STR", B4DIR, "/home/claude/torsion_jac/code");
RESDIR = if(type(RESDIR)=="t_STR", RESDIR, "/home/claude/torsion_jac/results");
read(strprintf("%s/claude_ov_b4_locus_data.gp", B4DIR));
if(type(WDATA)=="t_STR", read(WDATA));

PMAX = if(type(PMAX)=="t_INT", PMAX, 600);
TAG  = if(type(TAG)=="t_STR", TAG, "small");
DOW  = if(type(DOW)=="t_INT", DOW, 0);

u = 'u;
dlcm(P) = { my(L=1); for(i=0, poldegree(P), L = lcm(L, denominator(polcoef(P,i,u)))); L; };

/* build integral (p0,p1,d,c) with  value = c*(p0(u)+p1(u)*yy)*d(u)  in the
   same square class as (C0NUM/C0DEN) + (C1NUM/C1DEN)*yy               */
mkloc(C0N, C1N, DEN) = {
  my(L, M, p0, p1, d);
  L = lcm(dlcm(C0N), dlcm(C1N));
  M = dlcm(DEN);
  p0 = L*C0N; p1 = L*C1N; d = M*DEN;
  [p0, p1, d, L*M];
};

LOCA = mkloc(A1_C0_NUM, A1_C1_NUM, A1_C0_DEN);
LOCD = mkloc(D1_C0_NUM, D1_C1_NUM, D1_C0_DEN);
print("LOCA c=", LOCA[4], " degs=", [poldegree(LOCA[1]),poldegree(LOCA[2]),poldegree(LOCA[3])]);
print("LOCD c=", LOCD[4], " degs=", [poldegree(LOCD[1]),poldegree(LOCD[2]),poldegree(LOCD[3])]);
{
if(DOW,
  LOCW = mkloc(W_C0_NUM, W_C1_NUM, W_C0_DEN);
  print("LOCW c=", LOCW[4], " degs=", [poldegree(LOCW[1]),poldegree(LOCW[2]),poldegree(LOCW[3])]);
);
}

/* exact rational value of the test function at (uu,yyv) */
locval(LOC, uu, yyv) = LOC[4] * (subst(LOC[1],u,uu) + subst(LOC[2],u,uu)*yyv) * subst(LOC[3],u,uu);

/* ---- self-test: reproduce the direct-walk verdicts on the known members ---- */
{
  my(E, G, ok=1);
  E = ellinit([0,3,0,2,1]); G=[0,1];
  for(n=-40,-1,
    my(P=ellmul(E,G,n));
    if(#P<2, next);
    for(sg=0,1,
      my(uu=-P[1], yyv=(-1)^sg*P[2], a, d);
      if(uu==2 || uu==1, next);
      a = locval(LOCA, uu, yyv); d = locval(LOCD, uu, yyv);
      if(a != 0 && issquare(a), print("SELFTEST A1 SQUARE at n=",n," sg=",sg));
      if(d != 0 && issquare(d), print("SELFTEST D1 SQUARE at n=",n," sg=",sg));
      if(a == 0, print("SELFTEST A1 DEGEN (pole/zero) at n=",n," sg=",sg," u=",uu));
    ));
  print("SELFTEST done (expect exactly: D1 square at n=-6 sg=1, and A1/D1 DEGEN at u=-2)");
}

/* ---- table generation ---- */
{
  my(fn, plist, np=0, badp=0);
  fn = strprintf("%s/claude_ov_b4_tab_%s.txt", RESDIR, TAG);
  write(fn, "# p N_p ; char r <-> n = r mod N_p ; '1'=allowed (chi=+1 or 0), '0'=killed");
  plist = if(type(PLIST)=="t_VEC", PLIST, select(q->q>=5 && q!=23, primes([5,PMAX])));
  foreach(plist, p,
    my(Ep, Gp, Np, sA, sD, sW, Pt, bad=0);
    if(p<5 || p==23, next);
    Ep = ellinit([0,3,0,2,1], p);
    if(Ep == [], badp++; next);
    Gp = [Mod(0,p), Mod(1,p)];
    Np = ellorder(Ep, Gp);
    if(Np == 0 || Np == 1, badp++; next);
    sA = vector(Np); sD = vector(Np); if(DOW, sW = vector(Np));
    Pt = [Mod(0,p),Mod(1,p)];
    for(r = 1, Np,
      my(uu, yyv, a, d);
      if(r == Np, sA[Np]=1; sD[Np]=1; if(DOW,sW[Np]=1); break);  /* r=Np  <->  n=0 mod Np: P = O, allow */
      uu = -Pt[1]; yyv = Pt[2];
      a = locval(LOCA, uu, yyv); d = locval(LOCD, uu, yyv);
      sA[r] = if(a == 0, 1, if(issquare(a), 1, 0));
      sD[r] = if(d == 0, 1, if(issquare(d), 1, 0));
      if(DOW, my(ww = locval(LOCW, uu, yyv)); sW[r] = if(ww==0,1,if(issquare(ww),1,0)));
      Pt = elladd(Ep, Pt, Gp);
    );
    /* store with index r = n mod Np, r in [0,Np-1]; sA[r] above used r=n, so shift */
    my(strA, strD, strW);
    strA = concat(vector(Np, i, Str(sA[if(i==1, Np, i-1)])));
    strD = concat(vector(Np, i, Str(sD[if(i==1, Np, i-1)])));
    write(fn, Str(p, " ", Np));
    write(fn, strA);
    write(fn, strD);
    if(DOW,
      strW = concat(vector(Np, i, Str(sW[if(i==1, Np, i-1)])));
      write(fn, strW));
    np++;
    if(np % 25 == 0, print("PROGRESS primes=", np, " last p=", p));
  );
  print("MKTABLES_DONE tag=", TAG, " primes=", np, " skipped=", badp, " file=", fn);
}
quit
