\\ fpfile.gp — fingerprint points listed in fpin.txt (lines: p q rn rd label)
\\ against the two known curve classes; prints class or NEWCLASS per point.
default(parisize, "512M");
PL = [31, 37, 41, 43, 47, 53, 59, 61];
Av = [1,1,1,2,2];
mkcurve(pp, qq, rn, rd) = {
  my(u = pp/qq, r = rn/rd, s = u*(r-1), m = 4*r*u*(u-1), n = r-1, d, g, si, mi, ni, B);
  d = lcm([denominator(s), denominator(m), denominator(n)]);
  si = s*d; mi = m*d; ni = n*d;
  g = gcd([si, mi, ni]); si /= g; mi /= g; ni /= g;
  B = [2*si^2-si*ni, 2*si^2+si*mi-2*si*ni-mi*ni, 2*si^2+si*mi-si*ni-mi*ni, -mi*ni, 4*si^2-4*si*ni-mi*ni];
  [si, mi, ni, B];
}
fp(pp, qq, rn, rd) = {
  my(t = mkcurve(pp, qq, rn, rd), B = t[4], f, v = vector(#PL));
  f = prod(i=1, 5, Av[i] + B[i]*'x);
  for(j=1, #PL, my(pr = PL[j], cnt = 0, fr = f * Mod(1, pr));
    if(poldisc(f) % pr == 0 || pollead(f) % pr == 0, v[j] = -1,
      for(xx=0, pr-1, my(e = lift(subst(fr, 'x, xx)));
        cnt += 1 + kronecker(e, pr));
      cnt += 1;
      v[j] = cnt));
  v;
}
mycmp(v, w) = { my(ok = 1); for(j=1, #v, if(v[j] != -1 && w[j] != -1 && v[j] != w[j], ok = 0)); ok; }
{
c1 = fp(-97, 48, -49, 240);
c2 = fp(-23, 75, -9025, 3519);
V = readvec("fpin.txt");
for(i=1, #V,
  my(L = V[i], v = fp(L[1], L[2], L[3], L[4]));
  my(cl = if(mycmp(v, c1), "CURVE1", if(mycmp(v, c2), "CURVE2", "NEWCLASS")));
  printf("%s %s: fp=%s\n", cl, Str(L), Str(v)));
}
quit
