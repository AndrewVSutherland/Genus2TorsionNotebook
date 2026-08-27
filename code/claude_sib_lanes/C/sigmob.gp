\\ sigmob.gp — identify sigma geometrically: at random (u,r), match the 6
\\ Weierstrass points of model at (u,r) with those at sigma(u,r)=((4u-3)/(4u-4),1-r)
\\ by a Moebius map; report the induced slot permutation (slot 6 = infinity).
default(parisize, "512M");
Av = [1,1,1,2,2];
BS(s,m,n) = [2*s^2-s*n, 2*s^2+s*m-2*s*n-m*n, 2*s^2+s*m-s*n-m*n, -m*n, 4*s^2-4*s*n-m*n];
roots6(u,r) = {
  my(Bv = BS(u*(r-1), 4*r*u*(u-1), r-1));
  concat(vector(5, i, -Av[i]/Bv[i]), [oo]);
}
\\ Moebius through 3 points: M(z) with M(p1)=q1, M(p2)=q2, M(p3)=q3
\\ use cross-ratio: (M(z),q1;q2,q3) = (z,p1;p2,p3). Handle oo symbolically via
\\ projective coords: represent points as [a,b] (z=a/b), oo=[1,0].
prj(z) = if(z == oo, [1,0], [z,1]);
\\ cross ratio of 4 projective points
crp(A,B,C,D) = {
  my(d1 = A[1]*C[2]-A[2]*C[1], d2 = B[1]*D[2]-B[2]*D[1],
     d3 = A[1]*D[2]-A[2]*D[1], d4 = B[1]*C[2]-B[2]*C[1]);
  if(d3*d4 == 0, oo, (d1*d2)/(d3*d4));
}
{
uu = 3; rr = 5/7;
R1 = roots6(uu, rr);
us = (4*uu-3)/(4*uu-4); rs = 1-rr;
R2 = roots6(us, rs);
print("R1 = ", R1);
print("R2 = ", R2);
found = 0;
\\ try all injections of (R1[1],R1[2],R1[3]) onto ordered triples of R2;
\\ then check where remaining points go.
for(i=1,6, for(j=1,6, if(j==i, next); for(k=1,6, if(k==i||k==j, next);
  \\ candidate: R1[1]->R2[i], R1[2]->R2[j], R1[3]->R2[k]
  my(perm = vector(6, t, 0)); perm[1]=i; perm[2]=j; perm[3]=k;
  my(ok = 1);
  for(t=4,6,
    my(cr1 = crp(prj(R1[t]),prj(R1[1]),prj(R1[2]),prj(R1[3])));
    my(tgt = 0);
    for(l=1,6,
      if(l==i||l==j||l==k, next);
      my(cr2 = crp(prj(R2[l]),prj(R2[i]),prj(R2[j]),prj(R2[k])));
      if(cr1 == cr2, tgt = l; break));
    if(tgt == 0, ok = 0; break);
    \\ ensure injectivity
    for(tt=1,t-1, if(perm[tt]==tgt, ok=0));
    if(!ok, break);
    perm[t] = tgt);
  if(ok,
    found++;
    printf("MATCH perm (slot t of model1 -> slot perm[t] of model2): %s\n", Str(perm)));
)));
if(!found, print("NO Moebius matching — sigma-image has different branch config!"));
}
quit
