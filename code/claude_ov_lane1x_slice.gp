\\ claude_ov_lane1_slice.gp -- Lane 1: reduce the u = c slice of the contact-7
\\ three-root surface R(s,t,u)=0 to a plane cubic in (p,q)=(s+t,st) and study
\\ the resulting fibration.  Special attention to c = -1/2.
u='u; s='s; t='t; p='p; q='q; c='c; m='m; w='w;
A(w0) = (2*w0^5+4*w0^4+6*w0^3+8*w0^2+10*w0+5)/(2*(w0+1)^2);
{
bST = (A(t)-A(s))/(s^2-t^2);
aST = A(s) - bST*(1-s^2);
R = divrem(numerator(A(u) - aST - bST*(1-u^2)), (u-s)*(u-t))[1];
R = R/content(R);

\\ ---------- 1. rewrite R(s,t,c) in (p,q) ----------
\\ R is symmetric in s,t; express via power sums.
\\ Substitute t = p - s and reduce modulo s^2 - p s + q  (so that s,t are the roots).
Rc = subst(R, u, c);
\\ reduce Rc(s, p-s) modulo s^2-p*s+q as a polynomial in s
Rst = subst(Rc, t, p-s);
Rst = Rst % (s^2 - p*s + q);     \\ polynomial in s of degree <= 1
if(poldegree(Rst,s) > 0, print("NOT symmetric?  s-degree = ", poldegree(Rst,s)));
Rpq = subst(Rst, s, 0);
Rpq = Rpq/content(Rpq);
print("== R in (p,q,c) ==");
print("deg_p=",poldegree(Rpq,p)," deg_q=",poldegree(Rpq,q)," deg_c=",poldegree(Rpq,c));
print("Rpq = ", Rpq);
\\ total degree in (p,q)
print("");
print("-- terms of Rpq grouped by total (p,q)-degree --");
for(d=0,4, my(T=0);
  for(i=0,3, for(j=0,3, if(i+j==d, T += polcoef(polcoef(Rpq,i,p),j,q)*p^i*q^j)));
  if(T!=0, print("  deg ",d,": ",T)));

\\ ---------- 2. the c = -1/2 fibre ----------
Rm = subst(Rpq, c, -1/2); Rm = Rm/content(Rm);
print("");
print("== c = -1/2 plane cubic in (p,q) ==");
print("  ", Rm, "  = 0");
print("  value at (0,0): ", subst(subst(Rm,p,0),q,0));
print("  lowest-degree part at origin: ", polcoef(polcoef(Rm,1,p),1,q),"*p*q  (nondegenerate quadratic => NODE)");
\\ parametrize by p = m q
Par = subst(Rm, p, m*q);
Par = Par/q^2;
print("  after p = m*q, dividing by q^2: ", Par);
qsol = -polcoef(Par,0,q)/polcoef(Par,1,q);
print("  => q(m) = ", qsol, " = ", factor(numerator(qsol)), " / ", factor(denominator(qsol)));
psol = m*qsol;
print("  => p(m) = ", psol);
\\ splitting condition p^2-4q = square
D = psol^2 - 4*qsol;
print("  p^2-4q = ", D);
Dn = numerator(D); Dd = denominator(D);
print("     numerator   = ", factor(Dn));
print("     denominator = ", factor(Dd));

\\ ---------- 3. the five known slice points ----------
print("");
print("== the five known u=-1/2 points, in m ==");
sl = [[-10,-10/7],[-15/8,-15/19],[-511/61,-511/625],[-164/297,164/361],[-13/49,13/50]];
for(k=1,#sl, my(S0=sl[k][1], T0=sl[k][2], P0, Q0, M0, W2);
  P0 = S0+T0; Q0 = S0*T0; M0 = P0/Q0;
  W2 = -M0^4+6*M0^2+4*M0;
  print("  (s,t)=",sl[k],"  p=",P0," q=",Q0,"  m=",M0,"  -m^4+6m^2+4m = ",W2,
        "  issquare: ", issquare(W2)));
print("");
print("== the elliptic quartic  w^2 = -m^4 + 6 m^2 + 4 m = -m(m+2)(m^2-2m-2) ==");
print("  check factorisation: ", factor(-m^4+6*m^2+4*m));
}
quit
