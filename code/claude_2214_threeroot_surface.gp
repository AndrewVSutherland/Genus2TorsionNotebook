u='u; s='s; t='t;   \\ u highest priority
A(w) = (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2);
{
bST = (A(t)-A(s))/(s^2-t^2);
aST = A(s) - bST*(1-s^2);
E = A(u) - aST - bST*(1-u^2);
N = numerator(E);
dv = divrem(N, (u-s)*(u-t));
if(dv[2] != 0, print("UNEXPECTED remainder"); quit);
R = dv[1];
R = R/content(R);
print("deg_u=", poldegree(R,u), " deg_s=", poldegree(R,s), " deg_t=", poldegree(R,t));
\\ symmetry checks
Rc = substvec(R, [s,t,u], [t,u,s]);
Rs = substvec(R, [s,t,u], [t,s,u]);
print("cyclic-sym: ", Rc == R || Rc == -R, "  transposition-sym: ", Rs == R || Rs == -R);
print("R = ", R);
\\ verify the four hit orbits
hits = [[-3,-3/4,-3/5],[-10,-1/2,-10/7],[-5,-15/8,-15/22],[-1/2,-15/8,-15/19]];
for(k=1,4, print("hit",k," R=", substvec(R, [s,t,u], hits[k])));
\\ boundary factors: R at u=-1, s=t, u=s
print("R(u=-1) factored: ", factor(subst(R,u,-1)));
print("R(u=1): ", factor(subst(R,u,1)));
print("R(u=0): ", factor(subst(R,u,0)));
}
u='u; s='s; t='t;
A(w) = (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2);
{
bST = (A(t)-A(s))/(s^2-t^2);
aST = A(s) - bST*(1-s^2);
N = numerator(A(u) - aST - bST*(1-u^2));
R = divrem(N, (u-s)*(u-t))[1];
R = R/content(R);
\\ ---- fit R as polynomial in e1,e2,e3 ----
\\ basis: monomials e1^a e2^b e3^c, small degrees
bas = List();
for(a=0,6, for(b=0,3, for(c=0,3,
  if(a+2*b+3*c <= 9 && b+c <= 3 && a+b+c <= 6, listput(bas,[a,b,c])))));
bas = Vec(bas);
\\ evaluate at random rational triples
M = matrix(0,#bas); rhs = [];
count = 0; trial = 0;
while(count < #bas + 12 && trial < 400,
  trial++;
  my(sv=random(60)-30, tv=random(60)-30, uv=random(60)-30);
  if(#Set([sv,tv,uv]) != 3, next);
  my(e1=sv+tv+uv, e2=sv*tv+sv*uv+tv*uv, e3=sv*tv*uv);
  my(row = vector(#bas, k, e1^bas[k][1] * e2^bas[k][2] * e3^bas[k][3]));
  M = matconcat([M; Mat(row)]);
  rhs = concat(rhs, [substvec(R,[s,t,u],[sv,tv,uv])]);
  count++);

K = matker(matconcat([M, Mat(rhs~)]), 1);
if(#K == 0, print("no e-basis relation found (unexpected)"),
  cv = K[,1];
  den = -cv[#bas+1];
  print("R in e-basis (coeff * e1^a e2^b e3^c):");
  for(k=1, #bas,
    if(cv[k] != 0, print("  ", cv[k]/den, " * e1^",bas[k][1]," e2^",bas[k][2]," e3^",bas[k][3]))));
}
quit
