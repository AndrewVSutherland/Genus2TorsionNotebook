x='x;
SqrtPolyPart(f) = {
 my(s = x^3);
 for(k=1, 3,
  my(d = f - s^2);
  if(poldegree(d) <= 2, break);
  s += polcoef(d, 6-k)/(2*polcoef(s,3)) * x^(3-k));
 s }
CFPattern(f, maxsteps) = {
 my(s = SqrtPolyPart(f), Pi = 0, Qi = 1, pat = []);
 for(i=0, maxsteps,
  if(Qi == 0, return([pat, -1]));
  my(ai = divrem(Pi+s, Qi)[1]);
  pat = concat(pat, poldegree(ai));
  my(Pn = ai*Qi - Pi);
  if((f - Pn^2) % Qi != 0, return([pat, -2]));
  my(Qn = divrem(f - Pn^2, Qi)[1]);
  Pi = Pn; Qi = Qn;
  if(i >= 1 && poldegree(Qi) <= 0 && Qi != 0, return([pat, vecsum(pat)])));
 [pat, 0] }
{
F0 = 4*(9*x^6+69*x^5+123*x^4-95*x^3-183*x^2+165*x-35) + (x^2+x)^2;
fa = factor(F0);
rts = []; quad = 0;
for(k=1, matsize(fa)[1],
  g = fa[k,1];
  if(poldegree(g)==1, rts = concat(rts, [-polcoef(g,0)/polcoef(g,1)]));
  if(poldegree(g)==2, quad = g));
print("roots: ", rts);
print("quad: ", quad);
for(i=1, #rts, for(j=1, #rts,
  if(i==j, next);
  r1 = rts[i]; r2 = rts[j];
  Ft = subst(F0, x, r1 + (r2-r1)*x);
  Ft = Ft/pollead(Ft);
  pr = CFPattern(Ft, 20);
  oth = [];
  for(k=1, #rts, if(k!=i && k!=j, oth = concat(oth, [(rts[k]-r1)/(r2-r1)])));
  qq = subst(quad, x, r1 + (r2-r1)*x); qq = qq/pollead(qq);
  print("(",r1,"->0, ",r2,"->1): a,b=",oth, " c,d=",polcoef(qq,1),",",polcoef(qq,0), "  pattern=", pr[1], " order=", pr[2]);
));
}
quit
