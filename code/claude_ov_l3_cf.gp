\\ ============================================================================
\\ Lane 3 (route B2'): [1,1,2,2] + ord(D_inf)=11  ->  generic [2,22]
\\ Core machinery: exact polynomial-CF order of D_inf, the BLP (c,d,r1,r2)
\\ parametrization, and the factored [1,1,2,2] chart.
\\ Self-test:  gp -q code/claude_ov_l3_cf.gp
\\ ============================================================================
x = 'x;

\\ polynomial part of sqrt(f), f monic sextic
sqrtpolypart(f) = {
  my(s = x^3, dd);
  for(k=1,3,
    dd = f - s^2;
    if(poldegree(dd) <= 2, break);
    s = s + (polcoeff(dd, 6-k)/(2*polcoeff(s,3)))*x^(3-k));
  s;
}

\\ exact order of D_inf = inf+ - inf- via polynomial continued fraction.
\\ returns 0 if no quasi-period found within maxsteps.
cforder(f, maxsteps) = {
  my(s, Pi, Qi, total, ai, Pn, Qn);
  s = sqrtpolypart(f);
  Pi = 0*x; Qi = 1 + 0*x; total = 0;
  for(i=0, maxsteps,
    if(Qi == 0, return(0));
    ai = (Pi + s) \ Qi;
    total += poldegree(ai);
    Pn = ai*Qi - Pi;
    if((f - Pn^2) % Qi != 0, return(0));
    Qn = (f - Pn^2) \ Qi;
    Pi = Pn; Qi = Qn;
    if(i >= 1 && poldegree(Qi) <= 0 && Qi != 0, return(total)));
  0;
}

\\ ---- BLP2009 ansatz -------------------------------------------------------
\\ R = x^3 - x^2 + a x + b,  S = x^2 + d,  f = R^2 - 4 c^2 S^2 = (R-2cS)(R+2cS)
blp_f(a,b,c,d) = (x^3-x^2+a*x+b)^2 - 4*c^2*(x^2+d)^2;

\\ ROUTE B2' parametrization: the [1,1,2,2] condition is LINEAR in (a,b).
\\   (R-2cS)(r1) = 0 ,  (R+2cS)(r2) = 0   ->  solve for (a,b).
\\ returns [a,b]; requires r1 != r2.
blp_ab_from_roots(c,d,r1,r2) = {
  my(a,b);
  \\ r1^3 - (1+2c) r1^2 + a r1 + (b - 2cd) = 0
  \\ r2^3 + (2c-1) r2^2 + a r2 + (b + 2cd) = 0
  a = ( -(r1^3 - r2^3) + (1+2*c)*r1^2 + (2*c-1)*r2^2 + 4*c*d ) / (r1 - r2);
  b = -r1^3 + (1+2*c)*r1^2 - a*r1 + 2*c*d;
  [a,b];
}

\\ ---- factored [1,1,2,2] chart (the search chart) ---------------------------
\\ f = x (x-w) (x^2 + a x + b) (x^2 + c x + d)   -- monic sextic, 2-rank 2
chart_f(w,a,b,c,d) = x*(x-w)*(x^2+a*x+b)*(x^2+c*x+d);

nratroots(g) = { my(fa=factor(g), s=0); for(j=1,matsize(fa)[1], if(poldegree(fa[j,1])==1, s+=fa[j,2])); s; }
factortype(g) = { my(fa=factor(g)); vecsort(vector(matsize(fa)[1], j, poldegree(fa[j,1]))); }

{
  print("=== validated CF test vectors (pell-cf-order skill) ===");
  f14 = (x^2+1)*(x^4+5*x^2+4*x+4);
  f18 = (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8);
  f7  = x^6+2*x^5-5*x^4-14*x^3-3*x^2+24*x+28;
  printf("  f14 -> %d  (expect 14)\n", cforder(f14,40));
  printf("  f18 -> %d  (expect 18)\n", cforder(f18,40));
  printf("  f28cautionary -> %d  (expect 7)\n", cforder(f7,40));

  print("\n=== anchor 1: BLP corrected row C4 (a,b,c,d)=(1159/81,-277/243,40/9,13/27) ===");
  a=1159/81; b=-277/243; c=40/9; d=13/27;
  f = blp_f(a,b,c,d);
  printf("  CF order        = %d  (expect 11)\n", cforder(f,40));
  printf("  factor type     = %s\n", factortype(f));
  \\ recover (r1,r2) as the rational roots of the two cubics
  A = x^3-x^2+a*x+b - 2*c*(x^2+d); B = x^3-x^2+a*x+b + 2*c*(x^2+d);
  rA = [t[1] | t <- factor(A)~, poldegree(t[1])==1];
  print("  roots of R-2cS  : ", polrootsreal(A));
  fa=factor(A); r1 = -polcoeff(fa[1,1],0)/polcoeff(fa[1,1],1);
  fb=factor(B); r2 = -polcoeff(fb[1,1],0)/polcoeff(fb[1,1],1);
  printf("  r1 = %s   r2 = %s\n", r1, r2);
  ab = blp_ab_from_roots(c,d,r1,r2);
  printf("  ROUNDTRIP  (a,b) from (c,d,r1,r2) = (%s, %s)\n", ab[1], ab[2]);
  printf("  matches printed (a,b)?  %d\n", (ab[1]==a) && (ab[2]==b));
  f2 = blp_f(ab[1],ab[2],c,d);
  printf("  CF order of rebuilt f = %d\n", cforder(f2,40));
  \\ integral model x -> x/9 : f(x/9)*9^6
  g = subst(f, x, x/9)*9^6;
  printf("  integral model  = %s\n", g);
  printf("  == BLP printed integral model?  %d\n", g == (x-9)*(x+21)*(x^2-80*x+439)*(x^2+50*x+109));

  print("\n=== anchor 2: 19044.h.2 in the monic [1,1,2,2] chart (p)=(-3,8,4,27) ===");
  f19 = x*(x-1)*(x^2-3*x+8)*(x^2+4*x+27);
  printf("  CF order = %d  (expect 11)\n", cforder(f19,40));
  printf("  factor type = %s\n", factortype(f19));

  print("\n=== anchor 1 in the factored chart y^2 = x(x-w)(x^2+ax+b)(x^2+cx+d) ===");
  \\ C4corr normalized: (p1,p2,p3,p4)=(-34/15,32/45,31/15,-2/9), scale x by 15
  fc = x*(x-15)*(x^2-34*x+160)*(x^2+31*x-50);
  printf("  (w,a,b,c,d) = (15,-34,160,31,-50)  CF order = %d\n", cforder(fc,40));
  printf("  factor type = %s\n", factortype(fc));
}
quit
