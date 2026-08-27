\\ sieve22b.gp — fixed pilot sieve: extra rational root of the residual cubic
\\ on the two parametrized [22] families (=> 2-rank 2 => [2,22] candidate).
default(parisizemax, 2*10^9);
x = 'x; s = 's; T = 'T;
cA(ss) = (ss-1)^2*x^3 + (2*ss^4-4*ss^3+5*ss^2-4*ss+1)*x^2 \
        + (ss^6-2*ss^5+3*ss^4-2*ss^3+2*ss^2-2*ss+1)*x - ss^2*(ss^2-ss+1)^2;
Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;
teps(ss,e) = (-ss^2*(ss^2+1)*(ss^4-ss^2+1) + 2*e*ss^5)/((ss^2-1)^2);
{
F0 = Fl(teps(2,-1));
print("selftest A: ", (F0 % (x-4)) == 0, " ", (divrem(F0, x-4)[1] % cA(2)) == 0);
\\ family B cubic from Q4 (param log coefficients)
q3c = (211/100*T^3 - 3219/50*T^2 + 16417/25*T - 59954/25)/(T^3 - 23*T^2 + 168*T - 396);
q2c = (26041/40000*T^6 - 157789/10000*T^5 - 415573/2000*T^4 + 2879713/250*T^3 - 81235477/500*T^2 + 642935011/625*T - 1595049991/625)/(T^6 - 46*T^5 + 865*T^4 - 8520*T^3 + 46440*T^2 - 133056*T + 156816);
q1c = (-547353/500000*T^8 + 6482433/62500*T^7 - 66834977/15625*T^6 + 1578736769/15625*T^5 - 4708213621/3125*T^4 + 228428287724/15625*T^3 - 1415899642432/15625*T^2 + 5144130164272/15625*T - 8409021298024/15625)/(T^8 - 58*T^7 + 1453*T^6 - 20556*T^5 + 179820*T^4 - 997056*T^3 + 3425328*T^2 - 6671808*T + 5645376);
q0c = (194481/1000000*T^8 - 362502/15625*T^7 + 75029733/62500*T^6 - 550430994/15625*T^5 + 4013421667/6250*T^4 - 116565611824/15625*T^3 + 844813137732/15625*T^2 - 3501471367072/15625*T + 6374019201124/15625)/(T^8 - 58*T^7 + 1453*T^6 - 20556*T^5 + 179820*T^4 - 997056*T^3 + 3425328*T^2 - 6671808*T + 5645376);
dq = lcm([denominator(q3c),denominator(q2c),denominator(q1c),denominator(q0c)]);
Q4i = dq*(x^4 + q3c*x^3 + q2c*x^2 + q1c*x + q0c);
fa = factor(Q4i);
cB = 0;
for(k=1, matsize(fa)[1],
  if(poldegree(fa[k,1], x) == 3, cB = fa[k,1]));
print("family B cubic extracted: ", cB != 0);
}
mkpol(cc) = sum(j=1, #cc, cc[j]*x^(j-1));
mksieve(coefpols, plist) = {
  my(tabs = vector(#plist), pvar = variable(coefpols[1]));
  for(i=1, #plist,
    my(p = plist[i], tab = vector(p));
    for(r=0, p-1,
      my(cc = vector(#coefpols, j, subst(coefpols[j], pvar, r)));
      my(pol = mkpol(cc), ok);
      if(poldegree(pol) <= 0 || type(pol) != "t_POL", ok = 1,
        ok = (#polrootsmod(pol, p) > 0));
      tab[r+1] = ok);
    tabs[i] = tab);
  tabs;
}
runsieve(coefpols, H, plist, tag) = {
  my(tabs = mksieve(coefpols, plist), pvar = variable(coefpols[1]));
  my(nsurv = 0, nhit = 0, tot = 0);
  for(d=1, H, for(n=-H, H,
    if((n == 0 && d > 1) || gcd(abs(n), d) != 1, next);
    tot++;
    my(pass = 1);
    for(i=1, #plist,
      my(p = plist[i], dd = d % p);
      if(dd == 0, next);
      my(r = ((n % p) * lift(Mod(dd, p)^(-1))) % p);
      if(r < 0, r += p);
      if(!tabs[i][r+1], pass = 0; break));
    if(!pass, next);
    nsurv++;
    my(par = n/d);
    my(cc = vector(#coefpols, j, subst(coefpols[j], pvar, par)));
    my(pol = mkpol(cc));
    if(type(pol) != "t_POL" || poldegree(pol) < 3, next);
    my(rts = nfroots(, pol));
    if(#rts > 0, nhit++; print("HIT ", tag, " par=", par, " roots=", rts));
  ));
  print(tag, ": tot=", tot, " survivors=", nsurv, " hits=", nhit);
}
{
plist = primes([200, 460]); plist = plist[1..min(30,#plist)];
H = 1500;
cA0 = -s^2*(s^2-s+1)^2;
cA1 = s^6-2*s^5+3*s^4-2*s^3+2*s^2-2*s+1;
cA2 = 2*s^4-4*s^3+5*s^2-4*s+1;
cA3 = (s-1)^2;
runsieve([cA0,cA1,cA2,cA3], H, plist, "famA-");
runsieve([subst(cA0,s,-s),subst(cA1,s,-s),subst(cA2,s,-s),subst(cA3,s,-s)], H, plist, "famA+");
if(cB != 0,
  ccB = vector(4, j, polcoef(cB, j-1, x));
  runsieve(ccB, H, plist, "famB"));
print("SIEVE_PILOT_DONE");
}
quit
