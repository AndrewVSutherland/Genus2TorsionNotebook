x = 'x; aa='aa; bb='bb; cc='cc; dd4='dd4;
E0 = eval(readstr("/tmp/claude-1000/-home-claude-torsion-jac/482027be-f387-47ef-8672-ee258bea84e5/scratchpad/E0.txt")[1]);
SqrtPolyPart(f) = {
 my(s = x^3);
 for(k=1,3, my(d=f-s^2); if(poldegree(d)<=2, break); s += polcoef(d,6-k)/(2*polcoef(s,3))*x^(3-k));
 s }
CFPat(f, ms) = {
 my(s = SqrtPolyPart(f), Pi=0, Qi=1, pat=[]);
 for(i=0, ms,
  if(Qi==0, return([pat,-1]));
  my(ai = divrem(Pi+s, Qi)[1]);
  pat = concat(pat, poldegree(ai));
  my(Pn = ai*Qi-Pi);
  if((f-Pn^2)%Qi != 0, return([pat,-2]));
  my(Qn = divrem(f-Pn^2,Qi)[1]);
  Pi=Pn; Qi=Qn;
  if(i>=1 && poldegree(Qi)<=0 && Qi!=0, return([pat, vecsum(pat)])));
 [pat,0] }
{
anchors = [[-3/5,16/15,-4/15,-52/75],[-5/3,-16/9,4/9,-52/27],[15/16,-9/16,-1/4,-39/64],[8/5,-1/15,-26/15,1/25],[5/8,-1/24,-13/12,1/64],[-15,-24,26,9],[8/3,25/9,-22/9,-13/27],[3/8,25/24,-11/12,-13/192],[9/25,24/25,-22/25,-39/625],[1/16,25/16,-7/4,9/64],[16,25,-28,36],[16/25,1/25,-28/25,36/625]];
H = 60; nhit = 0; seen = Map();
sweep = List();
for(dn=1, H, for(nm=-H, H, if(gcd(abs(nm),dn)==1 && abs(nm/dn)<=H, listput(sweep, nm/dn))));
sweep = Vec(sweep);
print("sweep size per slice: ", #sweep);
for(ai=1, #anchors,
 A = anchors[ai];
 for(mode=1, 3,
  for(k=1, #sweep,
   v = sweep[k];
   my(a0,b0,c0);
   if(mode==1, a0=A[1]; b0=A[2]; c0=v);
   if(mode==2, a0=A[1]; b0=v;   c0=A[3]);
   if(mode==3, a0=v;   b0=A[2]; c0=A[3]);
   my(q5 = substvec(E0, [aa,bb,cc], [a0,b0,c0]));
   if(type(q5) != "t_POL" || poldegree(q5, dd4) < 1, next);
   my(rts = nfroots(, subst(q5, dd4, x)));
   foreach(rts, d0,
    my(key = Str([a0,b0,c0,d0]));
    if(mapisdefined(seen, key), next);
    mapput(seen, key, 1);
    if([a0,b0,c0,d0] == A, next);
    my(f = x*(x-1)*(x-a0)*(x-b0)*(x^2+c0*x+d0));
    if(poldegree(f) != 6 || poldisc(f) == 0, next);
    my(pr = CFPat(f, 20));
    if(pr[2] == 7,
     nhit++;
     print("S7HIT anchor#", ai, " mode", mode, " (a,b,c,d)=(", a0, ",", b0, ",", c0, ",", d0, ") pattern=", pr[1])));
  )));
print("total new order-7 points found: ", nhit);
}
quit
