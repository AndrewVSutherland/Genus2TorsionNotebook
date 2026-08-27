\\ probeA.gp - [2,2,14] via THIRD rational Weierstrass root on the contact-7
\\ two-root chart (see code/contact7_two_root_surface.m).  Roots r=1-v^2 of the
\\ quintic satisfy a + b*(1-v^2) = A(v)  (LINEAR in (a,b)),
\\ A(v) = (2v^5+4v^4+6v^3+8v^2+10v+5)/(2(v+1)^2).
\\ Sweep (s,t) -> (a,b); third root u <=> rational root of the residual cubic.
\\ Type [1,1,1,2] quintic => 2-rank 3 => torsion >= (Z/2)^3 x Z/7 = [2,2,14].
default(parisizemax, 2*10^9);
x = 'x; u = 'u;
H = 16;
Afun(w) = (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2);
goodp = [11,13,17,19,23,29,31,37,41,43,53,59];

{
rats = List();
for(d=1, H, for(n=-H, H, if(n!=0 && gcd(n,d)==1, listput(rats, n/d))));
rats = Vec(rats);
print("sweep values: ", #rats, "  pairs ~ ", #rats*(#rats-1)/2);
np = 0; nh3 = 0; ncand = 0;
for(i=1, #rats,
 s = rats[i];
 if(s == -1, next);
 for(j=i+1, #rats,
  t = rats[j];
  if(t == -1 || t == -s, next);
  np++;
  b = (Afun(t)-Afun(s))/(s^2-t^2);
  a = Afun(s) - b*(1-s^2);
  if(a+b == 5/2, next);
  N = 2*u^5+4*u^4+6*u^3+8*u^2+10*u+5 - (a + b*(1-u^2))*2*(u+1)^2;
  dv = divrem(N, (u-s)*(u-t));
  if(dv[2] != 0, print("DIVFAIL s=",s," t=",t); next);
  cub = dv[1];
  if(poldegree(cub) < 1, next);
  rts = nfroots(, cub);
  for(ri=1, #rts,
   uu = rts[ri];
   if(uu==s || uu==-s || uu==t || uu==-t || uu==-1 || uu==0, next);
   h = 1 - 7/2*x + a*x^2 + b*x^3;
   F7 = h^2 + (x-1)^7;
   dq = divrem(F7, x^2);
   if(dq[2] != 0, next);
   f = dq[1];
   if(poldegree(f) != 5 || poldisc(f) == 0, next);
   r1 = 1-s^2; r2 = 1-t^2; r3 = 1-uu^2;
   if(#Set([r1,r2,r3]) != 3, next);
   if(subst(f,x,r1) != 0 || subst(f,x,r2) != 0 || subst(f,x,r3) != 0,
      print("ROOTFAIL s=",s," t=",t," u=",uu); next);
   nh3++;
   fa = factor(f);
   degs = vecsort(apply(poldegree, fa[,1]~));
   gate56 = 1; gate112 = 1; gate168 = 1; rootpow = 0; ngood = 0;
   for(pi=1, #goodp,
    p = goodp[pi];
    chi = iferr(hyperellcharpoly(Mod(1,p)*f), E, 0);
    if(chi == 0 || poldegree(chi) != 4, next);
    ngood++;
    nJ = subst(chi, variable(chi), 1);
    if(nJ % 56 != 0, gate56 = 0);
    if(nJ % 112 != 0, gate112 = 0);
    if(nJ % 168 != 0, gate168 = 0);
    if(!rootpow && polisirreducible(chi),
     ok = 1;
     for(k=2, 12,
      if(!polisirreducible(charpoly(Mod(x^k, chi))), ok = 0; break));
     if(ok, rootpow = p)));
   print("HIT3 s=",s," t=",t," u=",uu," type=",degs,
         " ngood=",ngood," g56=",gate56," g112=",gate112,
         " g168=",gate168," rootpow_p=",rootpow);
   print("   a=",a," b=",b);
   print("   f = ", f);
   if(gate56 && rootpow, ncand++);
  );
 );
 if(i % 25 == 0, print("PROGRESS i=",i,"/",#rats," pairs=",np," hits3=",nh3," cands=",ncand));
);
print("DONE pairs=", np, " hits3=", nh3, " strong_candidates=", ncand);
}
quit
