\\ (4,3)-component point sweeps for Flynn and DS quadratic-factor loci,
\\ plus DS (6,2) exact disc; members get the quartic-split test ([2,22] probe).
x='x; t='t; u='u; v='v;
Fl = x^6+2*x^5+(2*t+3)*x^4+2*x^3+(t^2+1)*x^2+2*t*(1-t)*x+t^2;
DSf = x^6-4*x^5+8*(1+t)*x^4-(10+32*t)*x^3+8*(1+6*t+2*t^2)*x^2-4*(1+6*t+16*t^2)*x+64*t^2+1;
g43F = u^4 + (2*v - 2)*u^3 + (v^2 - 6*v + 3)*u^2 + (-8*v^2 + 10*v - 2)*u + (-4*v^3 + 8*v^2 - 4*v + 1);
g43D = 4*u^4 + (4*v + 8)*u^3 + (v^2 - 12*v + 8)*u^2 + (-16*v^2 - 30*v + 12)*u + (-4*v^3 - 16*v^2 - 12*v + 9);
{
\\ DS (6,2) exact disc
dv = divrem(DSf, x^2+u*x+v); r = dv[2];
e1 = polcoef(r,1,x); e0 = polcoef(r,0,x);
W = polresultant(e1,e0,t); fw = factor(W);
for(k=1, matsize(fw)[1],
  g = fw[k,1];
  if(poldegree(g,v)==2,
    dsc = polcoef(g,1,v)^2 - 4*polcoef(g,2,v)*polcoef(g,0,v);
    print("DS62 dsc EXACT = ", dsc)));
\\ (4,3) sweeps: for u rational, solve cubic-in-v; recover t; quartic-split test
HH = 48;
fams = [[Fl, g43F, "Flynn43"], [DSf, g43D, "DS43"]];
for(fi=1, 2,
 ff = fams[fi][1]; g43 = fams[fi][2]; nm = fams[fi][3];
 dvv = divrem(ff, x^2+u*x+v); rr = dvv[2];
 E1 = polcoef(rr,1,x); E0 = polcoef(rr,0,x);
 npt = 0;
 for(dd=1, HH, for(nn=-HH*dd, HH*dd,
  if(gcd(nn,dd) != 1, next);
  u0 = nn/dd;
  if(abs(u0) > HH, next);
  cubv = subst(g43, u, u0);
  if(poldegree(cubv, v) < 1, next);
  foreach(nfroots(, cubv), v0,
   if(v0 == 0, next);
   p1 = subst(subst(E1,u,u0),v,v0); p0 = subst(subst(E0,u,u0),v,v0);
   if(p1 == 0 && p0 == 0, next);
   g = if(p1==0, p0, if(p0==0, p1, gcd(p1,p0)));
   if(type(g) != "t_POL" || poldegree(g,t) < 1, next);
   foreach(nfroots(,g), t0,
    F = subst(ff, t, t0);
    if(poldegree(F) < 5 || poldisc(F) == 0, next);
    dq = divrem(F, x^2+u0*x+v0);
    if(dq[2] != 0, next);
    q4 = dq[1];
    fa = factor(q4);
    degs = vecsort(apply(poldegree, fa[,1]~));
    npt++;
    tag = if(degs==[2,2], "SPLIT22", if(degs==[1,1,2], "SPLIT112", if(degs==[1,1,1,1], "SPLITALL", "")));
    if(npt <= 60 || tag != "",
      print(nm, " u=",u0," v=",v0," t=",t0," q4=",degs," ",tag));
    if(tag != "", print("   F = ", F));
   );
  );
 ));
 print(nm, " members found: ", npt);
);
}
quit
