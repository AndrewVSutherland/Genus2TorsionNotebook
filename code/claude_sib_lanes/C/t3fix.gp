\\ t3fix.gp — tau on the T3 pencil ((12)-relabeling), fixed locus, and the
\\ condition collapse there.  T3 pencil: (s,m,n)=(uD,-uN,D), D=u-1+rho,
\\ N=(2-4rho)u-(1-4rho).  Conditions (mod squares, proven in prod_02):
\\ V1 == q(1-(4u-2)rho), V2 == (2u-1)q((2u-1)-(4u-4)rho), V4 == (u-1)q(rho+u-1),
\\ V5 == V1*V2*V4*(1-2rho), q = 4u^2-6u+3.
default(parisize, "2G");
samecl(f, g) = issquare(f*g);
sqcore(f) = {
  my(fa = factor(f), od = 1, cc);
  for(j=1, matsize(fa)[1], if(fa[j,2] % 2, od *= fa[j,1]));
  cc = simplify(f / factorback(fa));
  core(numerator(cc)*denominator(cc)) * od;
}
b = varhigher("b");
Av = [1,1,1,2,2];
BS(s,m,n) = [2*s^2-s*n, 2*s^2+s*m-2*s*n-m*n, 2*s^2+s*m-s*n-m*n, -m*n, 4*s^2-4*s*n-m*n];
{
DD = 'u - 1 + 'r;  NN = (2-4*'r)*'u - (1-4*'r);
Bv = BS('u*DD, -'u*NN, DD);
\\ tau = perm [2,1,3] fixing 4,5, affine b:
sg = [2,1,3,4,5];
BP = vector(5, i, Av[i] * Bv[sg[i]] * prod(j=1,5, if(j==sg[i], 1, b*Bv[j] - Av[j])));
lin = 2*BP[1]+2*BP[2]-2*BP[3]+BP[4]-BP[5];
quad = (BP[1]+BP[3]-BP[2])*BP[4] + 2*(BP[3]-BP[1]-BP[4])*(BP[3]-BP[2]);
g = gcd(lin, quad);
printf("T3 (12)-relabel: deg_b gcd = %d\n", poldegree(g, b));
if(poldegree(g, b) > 0,
  my(fa = factor(g));
  for(jf=1, matsize(fa)[1],
    my(ff = fa[jf,1]);
    if(poldegree(ff, b) == 1,
      my(b0 = -polcoef(ff, 0, b)/polcoef(ff, 1, b));
      my(BQ = vector(5, i, subst(BP[i], b, b0)));
      my(rp = BQ[3]/BQ[5]);   \\ T3 pencil coordinate rho = B3/B5
      my(up = ((BQ[1]+BQ[3]-BQ[2])/2)/(BQ[3]-BQ[2]));
      printf("  b0 = %s ; u' = %s ; rho' = %s\n", Str(b0), Str(simplify(up)), Str(simplify(rp))))));
\\ now the fixed locus of the found map and condition collapse:
\\ (fill in manually below once map known — compute generically)
q = 4*'u^2-6*'u+3;
V1 = q*(1-(4*'u-2)*'r);
V2 = (2*'u-1)*q*((2*'u-1)-(4*'u-4)*'r);
V4 = ('u-1)*q*('r+'u-1);
print("V1*V2 core: ", Str(sqcore(V1*V2)));
print("V1*V2*V4*(1-2r) core (V5 class): ", Str(sqcore(V1*V2*V4*(1-2*'r))));
}
quit
