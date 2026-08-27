\\ surf_id.gp — Lane A step 1: Pythagorean (u,g) chart of the T5 hit surface S.
\\ S: z^2 = (q*r-1)(q*r-2u+1), q = 4u^2-6u+3.  Writing h = q*r-u, w = u-1:
\\ z^2 + w^2 = h^2, so h = w(g^2+1)/(2g), z = w(1-g^2)/(2g) for rational g,
\\ i.e. r = r(u,g) = (2gu + (u-1)(g^2+1)) / (2g*q).
\\ Tasks: (a) verify X1*X2*X3*X4 is a perfect square in Q(u,g) [V'4 dependent];
\\ (b) square-class normal forms W1,W2,W3 (+W4) on the chart, with degrees;
\\ (c) validate W-classes against X-classes at random points;
\\ (d) recover g for the two known hit representations; show the 6 near-misses
\\     lie OFF S (z^2 nonsquare);  (e) emit integer-cleared forms for the C sweep.
default(parisize, "2G");
q(u) = 4*u^2 - 6*u + 3;
X1f(u,r) = r*(q(u)*r - (2*u-1));
X2f(u,r) = q(u)*r^2 - (4*u^2-4*u+2)*r + (2*u-1);
c3f(u) = 16*u^4-40*u^3+40*u^2-18*u+3;
c2f(u) = -16*u^4+32*u^3-28*u^2+10*u-1;
c1f(u) = 8*u^3-12*u^2+10*u-3;
c0f(u) = -2*u+1;
X3f(u,r) = r*(c3f(u)*r^3 + c2f(u)*r^2 + c1f(u)*r + c0f(u));
X4f(u,r) = c3f(u)*r^2 + (-16*u^3+28*u^2-18*u+4)*r + (2*u-1)^2;
sqclv(x) = core(numerator(x)*denominator(x));  \\ square class of a rational
\\ square class of a rational function in (u,g): signed squarefree kernel poly
oddpart(P) = {
  my(fa = factor(P), od = 1, cc);
  cc = simplify(P / prod(j=1, matsize(fa)[1], fa[j,1]^fa[j,2]));  \\ constant content
  if(type(cc) == "t_POL", cc = simplify(subst(subst(cc, g, 0), u, 0)));
  for(j=1, matsize(fa)[1], if(fa[j,2]%2, od *= fa[j,1]));
  core(cc)*od;
}
{
g = varhigher("g", u);
rug = (2*g*u + (u-1)*(g^2+1)) / (2*g*q(u));
\\ sanity: is (q*r-1)(q*r-2u+1) a square in Q(u,g)?
zz = (q(u)*rug-1)*(q(u)*rug-2*u+1);
print("S-membership check: (qr-1)(qr-2u+1) square in Q(u,g)? ", issquare(zz));
X = vector(4);
X[1] = X1f(u,rug); X[2] = X2f(u,rug); X[3] = X3f(u,rug); X[4] = X4f(u,rug);
print("X1*X2*X3*X4 square in Q(u,g)? ", issquare(X[1]*X[2]*X[3]*X[4]));
W = vector(4);
for(k=1,4,
  my(P = numerator(X[k])*denominator(X[k]));
  W[k] = oddpart(P);
  printf("W%d: deg_u=%d deg_g=%d\n  = %s\n  factor: %s\n",
    k, poldegree(W[k], u), poldegree(W[k], g), Str(W[k]), Str(factor(W[k]))));
print("W4 == W1*W2*W3 mod squares? ", issquare(W[1]*W[2]*W[3]*W[4]));
\\ (c) random validation: class(Xk(u0,r0)) == class(Wk(u0,g0)) at 60 points
my(bad = 0, cnt = 0);
for(t=1, 200,
  if(cnt >= 60, break);
  my(u0 = random(400)/(random(200)+1) - 1, g0 = (random(400)+1)/(random(200)+1) - 1);
  if(g0 == 0 || u0 == 0 || u0 == 1 || 2*u0 == 1 || g0 == 1 || g0 == -1, next);
  my(r0 = subst(subst(rug, g, g0), u, u0));
  if(r0 == 0 || r0 == 1, next);
  my(ok = 1);
  for(k=1,4,
    my(Xv = if(k==1,X1f(u0,r0),if(k==2,X2f(u0,r0),if(k==3,X3f(u0,r0),X4f(u0,r0)))));
    my(Wv = subst(subst(W[k], g, g0), u, u0));
    if(Xv == 0 || Wv == 0, ok = -1; break);
    if(sqclv(Xv) != sqclv(Wv), ok = 0; break));
  if(ok == 0, bad++; printf("  MISMATCH at u=%s g=%s\n", Str(u0), Str(g0)));
  if(ok == 1, cnt++));
printf("random validation: %d points ok, %d mismatches\n", cnt, bad);
\\ (d) recover g at the two hit representations; near-miss off-S proof
hits = [[-97/48, -49/240], [133/145, 289/240]];
for(i=1,2,
  my(u0 = hits[i][1], r0 = hits[i][2], h = q(u0)*r0 - u0, w = u0 - 1);
  my(z2 = h^2 - w^2);
  my(z); if(!issquare(z2, &z), printf("hit %d NOT on S?!\n", i),
    my(g1 = (h+z)/w, g2 = (h-z)/w);
    printf("hit %d: u=%s r=%s  g = %s or %s (product %s)\n",
      i, Str(u0), Str(r0), Str(g1), Str(g2), Str(g1*g2));
    \\ confirm r(u,g) round-trips
    for(j=1,2, my(gg = if(j==1,g1,g2));
      if(subst(subst(rug, g, gg), u, u0) != r0, printf("  ROUNDTRIP FAIL g=%s\n", Str(gg))))));
nears = [[-4, -169/1431], [10, -25/551], [17, -1/143], [13/4, -1/143], [41/25, -25/551], [43/52, 841/697]];
for(i=1,#nears,
  my(u0 = nears[i][1], r0 = nears[i][2], h = q(u0)*r0 - u0, w = u0-1);
  printf("near %d: u=%s r=%s  z^2 class = %s (on S iff 1)\n",
    i, Str(u0), Str(r0), Str(sqclv(h^2 - w^2))));
}
quit
