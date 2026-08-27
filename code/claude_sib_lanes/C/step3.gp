\\ step3.gp — group structure of <sigma,tau>, fix(tau) condition collapse,
\\ (45)-composed chart matchings (T4 subsumption test).
default(parisize, "2G");
samecl(f, g) = issquare(f*g);
sqcore(f) = {
  my(fa = factor(f), od = 1, cc);
  for(j=1, matsize(fa)[1], if(fa[j,2] % 2, od *= fa[j,1]));
  cc = simplify(f / factorback(fa));
  core(numerator(cc)*denominator(cc)) * od;
}
\\ maps in (u,r) coordinates
sigU(u,r) = [(4*u-3)/(4*u-4), 1-r];
tauU(u,r) = [(4*r*u-3*r-1)/(4*r*(u-1)), r];
b = varhigher("b");
Av = [1,1,1,2,2];
BS(s,m,n) = [2*s^2-s*n, 2*s^2+s*m-2*s*n-m*n, 2*s^2+s*m-s*n-m*n, -m*n, 4*s^2-4*s*n-m*n];
{
\\ ---- group structure ----
st = sigU(tauU('u,'r)[1], tauU('u,'r)[2]);
ts = tauU(sigU('u,'r)[1], sigU('u,'r)[2]);
printf("sigma tau = [%s, %s]\n", Str(st[1]), Str(st[2]));
printf("tau sigma = [%s, %s]\n", Str(ts[1]), Str(ts[2]));
stst = sigU(tauU(st[1], st[2])[1], tauU(st[1], st[2])[2]);
printf("(sigma tau)^2 = [%s, %s]\n", Str(simplify(stst[1])), Str(simplify(stst[2])));
\\ ---- fix(tau): r = -1/(v(v-2)), v = 2u-1 -> 4u^2-8u+3 = (2u-1)(2u-3) ----
print("=== fix(tau) condition collapse (verify) ===");
Q = 'v^2-'v+1;
rfix = -1/('v*('v-2));
X = Q*rfix;
F = [Q*X*(X-'v), Q*(X-Q)*(X-'v), 'v*Q*X*(X-1)*(X-Q)*('v*X-Q), 'v*Q*(X-'v)*('v*X-Q)];
printf("V'1 == v+1 on fix: %d\n", samecl(F[1], 'v+1));
printf("V'2 == v+1 on fix: %d\n", samecl(F[2], 'v+1));
printf("V'3 == (2v-1)Q on fix: %d\n", samecl(F[3], (2*'v-1)*Q));
printf("V'4 == v^2-1 on fix: %d\n", samecl(F[4], ('v-1)*('v+1)));
\\ with v = (w^4+4)/(4w^2) (making v+1 and v-1 squares):
vw = ('w^4+4)/(4*'w^2);
printf("v+1 at vw is square: %d\n", issquare(vw+1));
printf("v-1 at vw is square: %d\n", issquare(vw-1));
c3w = subst((2*'v-1)*Q, 'v, vw);
printf("remaining condition core: %s\n", Str(sqcore(c3w)));
printf("factored: %s\n", Str(factor(sqcore(c3w))));
\\ ---- (45)-composed relabelings: does T4 pull back to T5? ----
print("=== chart matchings with slots 4,5 swapped ===");
Bv = BS('u*('r-1), 4*'r*'u*('u-1), 'r-1);
perms = [[1,2,3],[2,1,3],[1,3,2],[3,2,1],[2,3,1],[3,1,2]];
for(ip=1, 6,
  my(sg = concat(perms[ip], [5,4]));
  my(BP = vector(5, i,
    Av[i] * Bv[sg[i]] * prod(j=1,5, if(j==sg[i], 1, b*Bv[j] - Av[j]))));
  my(lin = 2*BP[1]+2*BP[2]-2*BP[3]+BP[4]-BP[5]);
  my(quad = (BP[1]+BP[3]-BP[2])*BP[4] + 2*(BP[3]-BP[1]-BP[4])*(BP[3]-BP[2]));
  my(g = gcd(lin, quad));
  printf("perm %s +swap45: deg_b gcd = %d\n", Str(perms[ip]), poldegree(g, b));
  if(poldegree(g, b) > 0,
    my(fa = factor(g));
    for(jf=1, matsize(fa)[1],
      my(ff = fa[jf,1]);
      if(poldegree(ff, b) == 1,
        my(b0 = -polcoef(ff, 0, b)/polcoef(ff, 1, b));
        my(BQ = vector(5, i, subst(BP[i], b, b0)));
        my(rp = BQ[4]/BQ[5]);
        my(up = ((BQ[1]+BQ[3]-BQ[2])/2)/(BQ[3]-BQ[2]));
        printf("  b0 = %s ; u' = %s ; r' = %s\n", Str(b0), Str(up), Str(rp))))));
}
quit
