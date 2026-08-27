x = 'x;
p = 101;
SqrtPolyPartP(f) = {
 my(s = Mod(1,p)*x^3);
 for(k=1,3, my(d=f-s^2); if(poldegree(d)<=2, break); s += polcoef(d,6-k)/(2*polcoef(s,3))*x^(3-k));
 s }
CFPatP(f, ms) = {
 my(s = SqrtPolyPartP(f), Pi=0, Qi=Mod(1,p), pat=[]);
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
cnt = Map(); tot = 0;
for(i=1, 400000,
 my(aa=random(p), bb=random(p), cc=random(p), dd=random(p));
 my(f = Mod(1,p)*x*(x-1)*(x-aa)*(x-bb)*(x^2+cc*x+dd));
 if(poldegree(f)!=6 || poldegree(gcd(f,f'))>0, next);
 my(pr = CFPatP(f, 20));
 if(pr[2] == 7,
  tot++;
  my(k = Str(pr[1]), c = if(mapisdefined(cnt,k), mapget(cnt,k), 0));
  mapput(cnt, k, c+1)));
print("order-7 hits: ", tot);
print("pattern census: ", Mat(cnt));
}
quit
