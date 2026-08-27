\\ Daowsud-Schmidt arXiv:1708.05511v2 (corrected) order-11 family: first probe.
x = 'x;
g(u) = x^6 - (16*u/(u^5+8))*x^5 \
 - ((u^15+40*u^10+128*u^5+512)/(2*u^3*(u^5+8)^2))*x^4 \
 + ((6*u^15+224*u^10+1536*u^5+3072)/(u^2*(u^5+8)^3))*x^3 \
 + ((u^30+112*u^25+2880*u^20+25600*u^15+106496*u^10+262144*u^5+262144)/(16*u^6*(u^5+8)^4))*x^2 \
 - ((u^25+80*u^20+1664*u^15+12288*u^10+36864*u^5+32768)/(2*u^5*(u^5+8)^4))*x \
 - (u^25+46*u^20+736*u^15+5248*u^10+18432*u^5+24576)/(2*u^4*(u^5+8)^4);
SqrtPolyPart(f) = {
 my(s = x^3);
 for(k=1, 3,
  my(d = f - s^2);
  if(poldegree(d) <= 2, break);
  s += polcoef(d, 6-k)/(2*polcoef(s,3)) * x^(3-k));
 s }
CFOrd(f, maxsteps) = {
 my(s = SqrtPolyPart(f), Pi = 0, Qi = 1, total = 0);
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
{
\\ transcription check: CF order must be 11
foreach([1, 2, 1/2, -3, 5, -1/4], u0,
  print("v2 u=", u0, " CFOrd=", CFOrd(g(u0), 60)));
\\ factor-type sweep
HH = 32; types = Map(); interesting = 0;
for(d=1, HH, for(n=-HH*d, HH*d,
 if(n == 0 || gcd(n,d) != 1 || abs(n/d) > HH, next);
 my(u0 = n/d);
 my(f = g(u0));
 if(type(f) != "t_POL" || poldegree(f) != 6 || poldisc(f) == 0, next);
 my(fa = factor(f), degs = vecsort(apply(poldegree, fa[,1]~)));
 my(k = Str(degs), c = if(mapisdefined(types,k), mapget(types,k), 0));
 mapput(types, k, c+1);
 if(degs != [6] && degs != [1,5],
   interesting++;
   if(interesting <= 25, print("V2MEMBER u=", u0, " type=", degs)));
));
print("type census: ", Mat(types));
}
quit
