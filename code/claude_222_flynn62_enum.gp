x='x; t='t; u='u; v='v;
Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;
\\ (6,2)-component data: a2*v^2 + b1*v + c0 with
a2(uu) = (uu-2)^4;
b1(uu) = 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);
c0(uu) = uu^2*(uu^2-uu+1)^2;
\\ E: Y^2 = w^3+3w^2+2w+1, u = -w; sqrt(dsc) = 4(u-1)*Y
{
E = ellinit([0,3,0,2,1]);
G = [0,1];  \\ candidate generator
pts = List();
P = G;
for(k=1, 14, listput(pts, P); listput(pts, ellneg(E,P)); P = elladd(E, P, G));
\\ also throw in small search points
foreach(hyperellratpoints(x^3+3*x^2+2*x+1, 3000), q, if(#q==2, listput(pts, q)));
seen = Map();
nmem = 0;
foreach(Vec(pts), Pt,
  w = Pt[1]; Y = Pt[2];
  uu = -w;
  if(uu == 2 || uu == 1, next);
  A = a2(uu); B = b1(uu);
  for(sg=0, 1,
    vv = (-B + (-1)^sg*4*(uu-1)*Y)/(2*A);
    if(vv == 0, next);
    \\ recover t: common root of the two remainder coefficients
    ff = Fl(t);
    dv = divrem(ff, x^2+uu*x+vv);
    r = dv[2];
    ee1 = polcoef(r,1,x); ee0 = polcoef(r,0,x);
    if(ee1 == 0 && ee0 == 0, next);
    g = if(ee1==0, ee0, if(ee0==0, ee1, gcd(ee1,ee0)));
    if(type(g) != "t_POL" || poldegree(g,t) < 1, next);
    foreach(nfroots(,g), t0,
      key = Str([uu,vv,t0]);
      if(mapisdefined(seen,key), next);
      mapput(seen,key,1);
      F = subst(ff, t, t0);
      if(poldegree(F) < 5 || poldisc(F) == 0, next);
      dq = divrem(F, x^2+uu*x+vv);
      if(dq[2] != 0, print("NODIV u=",uu," v=",vv," t=",t0); next);
      q4 = dq[1];
      fa = factor(q4);
      degs = vecsort(apply(poldegree, fa[,1]~));
      nmem++;
      splitflag = if(degs == [2,2], "SPLIT22", if(degs==[1,1,2], "SPLIT112", if(degs==[1,1,1,1],"SPLITALL", if(degs==[1,3],"SPLIT13",""))));
      print("MEMBER u=",uu," v=",vv," t=",t0," q4type=",degs," ",splitflag);
      if(splitflag != "", print("   F = ", F));
    );
  );
);
print("total distinct members: ", nmem);
}
quit
