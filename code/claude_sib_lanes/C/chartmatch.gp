\\ chartmatch.gp — find ALL model self-matchings of the T5 pencil under
\\ relabelings sigma of the three A=1 slots {1,2,3} (class-T5-preserving),
\\ with affine Moebius x -> alpha x + beta (slot 6 = infinity fixed).
\\ Chart: A=[1,1,1,2,2], B_i(s,m,n) quadrics; relations:
\\   (lin)  B5 = 2B1+2B2-2B3+B4
\\   (quad) (B1+B3-B2)*B4 + 2*(B3-B1-B4)*(B3-B2) = 0   [(s^2)(mn)=(sm)(sn)]
\\ Recovery: (s:m:n) = ((B1+B3-B2)/2 : B3-B1-B4 : B3-B2).
\\ On the pencil: s=u(r-1), m=4ru(u-1), n=r-1;  u'=s'/n', r'=B4'/B5'.
\\ New model roots: -A_i/B'_i = x_{sigma(i)} + b  (b = beta/alpha, one param).
\\ B'_i \propto A_i * B_{sigma(i)} * prod_{j != sigma(i)} (b*B_j - A_j).
default(parisize, "2G");
b = varhigher("b");
Av = [1,1,1,2,2];
BS(s,m,n) = [2*s^2-s*n, 2*s^2+s*m-2*s*n-m*n, 2*s^2+s*m-s*n-m*n, -m*n, 4*s^2-4*s*n-m*n];
{
Bv = BS('u*('r-1), 4*'r*'u*('u-1), 'r-1);
perms = [[1,2,3],[2,1,3],[1,3,2],[3,2,1],[2,3,1],[3,1,2]];
for(ip=1, 6,
  my(sg = concat(perms[ip], [4,5]));
  my(BP = vector(5, i,
    Av[i] * Bv[sg[i]] * prod(j=1,5, if(j==sg[i], 1, b*Bv[j] - Av[j]))));
  my(lin = 2*BP[1]+2*BP[2]-2*BP[3]+BP[4]-BP[5]);
  my(quad = (BP[1]+BP[3]-BP[2])*BP[4] + 2*(BP[3]-BP[1]-BP[4])*(BP[3]-BP[2]));
  my(g = gcd(lin, quad));
  printf("perm %s: deg_b gcd = %d\n", Str(perms[ip]), poldegree(g, b));
  if(poldegree(g, b) > 0,
    my(fa = factor(g));
    for(jf=1, matsize(fa)[1],
      my(ff = fa[jf,1]);
      if(poldegree(ff, b) == 1,
        my(b0 = -polcoef(ff, 0, b)/polcoef(ff, 1, b));
        \\ evaluate BP at b0
        my(BQ = vector(5, i, subst(BP[i], b, b0)));
        my(rp = BQ[4]/BQ[5]);
        my(up = ((BQ[1]+BQ[3]-BQ[2])/2)/(BQ[3]-BQ[2]));
        printf("  b0 = %s\n", Str(b0));
        printf("  u' = %s\n", Str(up));
        printf("  r' = %s\n", Str(rp));
        \\ numeric check at (u,r) = (17, -1/143):
        my(u0 = 17, r0 = -1/143);
        my(upn = substvec(up, ['u,'r], [u0,r0]), rpn = substvec(rp, ['u,'r], [u0,r0]));
        printf("  at (17,-1/143) -> (%s, %s)\n", Str(upn), Str(rpn));
      ,
        printf("  (degree-%d factor in b: %s)\n", poldegree(ff,b), Str(ff))))));
}
quit
