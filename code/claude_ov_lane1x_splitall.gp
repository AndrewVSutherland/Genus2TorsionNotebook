\\ claude_ov_lane1_splitall.gp -- Lane 1: the ORDER-112 ([2,2,2,14]) condition
\\ restricted to the u = -1/2 family.
\\ Given three roots v1=-1/2, v2=s, v3=t of the monic quintic
\\    Q(v) = v^5 + c4 v^4 + (2c4-1)v^3 + (c4+c0-1/2)v^2 + 2c0 v + c0,
\\ the residual quadratic is v^2 + Bv + C with
\\    B = c4 + (s+t+u),   C = -c0/(s t u),   D2 = B^2-4C.
\\ SPLITALL (5 rational roots, 2-rank 4, torsion [2,2,2,14]) <=> D2 is a square.
\\ On the family  p=s+t=2m^2/D, q=st=2m/D, D=(m+1)^2(m-2),  D2 is a rational
\\ function of m alone.  Compute it.
m='m; p='p; q='q; s='s; t='t;
G(v) = -(v^5 - v^3 - v^2/2)/(v+1)^2;
{
\\ ---- c4, c0 as symmetric functions of (s,t) ----
c4 = (G(s)-G(t))/(s^2-t^2);
c0 = G(s) - c4*s^2;
\\ verify against the direct chart on a numeric point
S0=-13/49; T0=13/50;
print("numeric check c4,c0 at (s,t)=(",S0,",",T0,"): ",
      subst(subst(c4,s,S0),t,T0), "  ", subst(subst(c0,s,S0),t,T0));

\\ ---- residual quadratic discriminant, symbolically in (s,t) ----
u0 = -1/2;
B = c4 + (s+t+u0);
C = -c0/(s*t*u0);
D2 = B^2 - 4*C;
\\ sanity at the numeric point: rebuild the quintic and divide out the 3 roots
c4n = subst(subst(c4,s,S0),t,T0); c0n = subst(subst(c0,s,S0),t,T0);
v='v;
Qn = v^5 + c4n*v^4 + (2*c4n-1)*v^3 + (c4n+c0n-1/2)*v^2 + 2*c0n*v + c0n;
resq = Qn/((v-S0)*(v-T0)*(v-u0));
print("residual quadratic at that point: ", resq);
print("  its discriminant: ", polcoef(resq,1)^2 - 4*polcoef(resq,0));
print("  D2 formula gives: ", subst(subst(D2,s,S0),t,T0));

\\ ---- push D2 into (p,q) then into m ----
\\ D2 is symmetric in s,t: reduce numerator/denominator mod s^2-p*s+q
D2n = numerator(D2); D2d = denominator(D2);
red(P0) = my(A); A = subst(P0, t, p-s) % (s^2-p*s+q); if(poldegree(A,s)>0, print("NOT SYMMETRIC")); subst(A,s,0);
D2npq = red(D2n); D2dpq = red(D2d);
D2pq = D2npq/D2dpq;
print("");
print("D2 in (p,q) = ", D2pq);

\\ ---- substitute the elliptic parametrisation ----
Dm = (m+1)^2*(m-2);
pm = 2*m^2/Dm; qm = 2*m/Dm;
D2m = subst(subst(D2pq, p, pm), q, qm);
D2m = D2m + 0;      \\ force simplification
NN = numerator(D2m); DD = denominator(D2m);
print("");
print("D2(m) numerator   = ", factor(NN));
print("D2(m) denominator = ", factor(DD));
\\ squarefree part of D2 as a rational function
sqfree(F) = my(fa, res); fa = factor(F); res = 1;
  for(i=1,#fa[,1], if(fa[i,2]%2==1, res *= fa[i,1])); res;
SF = sqfree(NN)*sqfree(DD);
print("");
print("=> D2(m) is a square  <=>  z^2 = ", SF, "   (squarefree part, up to the constant)");
print("   degree = ", poldegree(SF));
print("   factored: ", factor(SF));
print("   leading coeff / content: ", content(SF), "  lc=", pollead(SF));

\\ ---- numeric confirmation on the known family members ----
print("");
print("== check D2(m) on the known family members ==");
ms = [-4/5,-9/5,-98/73,16/41,1/13];
for(i=1,#ms, my(M=ms[i], val, sf);
  val = subst(D2m, m, M);
  sf = subst(SF, m, M);
  print("  m=",M,"  D2=",val,"  issquare(D2)=", issquare(val), "  sqfreepart=",sf, "  issquare=", issquare(sf)));
}
quit
