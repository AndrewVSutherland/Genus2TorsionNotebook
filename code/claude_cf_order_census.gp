\\ probeC.gp - mod-p density census of CF (D_infinity) order loci:
\\  shape A: (x^2+b1x+c1)(x^2+b2x+c2)(x^2+b3x+c3)  -> orders 11, 22 ([2,22] side)
\\  shape B: x(x-1)(x-a)(x-b)(x^2+cx+d)            -> orders 7, 14 ([2,2,14] side)
\\ Codim-2 locus in the shape space => expect counts ~ c * N / p^2; measures c
\\ and confirms nonemptiness over many F_p (prereq for backward-CF routes).
default(parisizemax, 1*10^9);
x = 'x;
NS = 200000;

SqrtPolyPartP(f) = {
 my(s = Mod(1, p)*x^3);
 for(k=1, 3,
  my(d = f - s^2);
  if(poldegree(d) <= 2, break);
  s += polcoef(d, 6-k)/(2*polcoef(s,3)) * x^(3-k));
 s }
CFOrdP(f, maxsteps) = {
 my(s = SqrtPolyPartP(f), Pi = 0, Qi = Mod(1,p), total = 0);
 for(i=0, maxsteps,
  if(Qi == 0, return(0));
  my(ai = divrem(Pi+s, Qi)[1]);
  total += poldegree(ai);
  my(Pn = ai*Qi - Pi);
  if((f - Pn^2) % Qi != 0, return(0));
  my(Qn = divrem(f - Pn^2, Qi)[1]);
  Pi = Pn; Qi = Qn;
  if(i >= 1 && poldegree(Qi) <= 0 && Qi != 0, return(total)));
 0 }

{for(pi=1, 2,
 p = [101, 211][pi];
 \\ shape A
 my(cnt = vector(60), good = 0, seeds11 = List());
 for(i=1, NS,
  my(b1=random(p), c1=random(p), b2=random(p), c2=random(p), b3=random(p), c3=random(p));
  my(f = Mod(1,p)*(x^2+b1*x+c1)*(x^2+b2*x+c2)*(x^2+b3*x+c3));
  if(poldegree(gcd(f, f')) > 0, next);
  good++;
  my(o = CFOrdP(f, 45));
  if(o >= 1 && o <= 59, cnt[o+1]++);
  if(o == 11 && #seeds11 < 4, listput(seeds11, [b1,c1,b2,c2,b3,c3])));
 print("p=", p, " shapeA[2,2,2] good=", good);
 foreach([7,11,13,14,22], o,
  print("   ord", o, ": ", cnt[o+1], "  c=", (cnt[o+1]*p^2*1.)/good));
 print("   seeds11: ", Vec(seeds11));
 \\ shape B
 my(cntB = vector(60), goodB = 0, seeds7 = List());
 for(i=1, NS,
  my(aa=random(p), bb=random(p), cc=random(p), dd=random(p));
  my(f = Mod(1,p)*x*(x-1)*(x-aa)*(x-bb)*(x^2+cc*x+dd));
  if(poldegree(f) != 6 || poldegree(gcd(f, f')) > 0, next);
  goodB++;
  my(o = CFOrdP(f, 45));
  if(o >= 1 && o <= 59, cntB[o+1]++);
  if(o == 7 && #seeds7 < 4, listput(seeds7, [aa,bb,cc,dd])));
 print("p=", p, " shapeB[1,1,1,1,2] good=", goodB);
 foreach([5,7,9,11,14], o,
  print("   ord", o, ": ", cntB[o+1], "  c=", (cntB[o+1]*p^2*1.)/goodB));
 print("   seeds7: ", Vec(seeds7));
);}
print("PROBE_C_DONE");
quit
