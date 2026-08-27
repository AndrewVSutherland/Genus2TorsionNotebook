\\ claude_ov_lane1_fibration.gp -- Lane 1: the u = c slices of R(s,t,u)=0.
\\ In (p,q)=(s+t,st) each slice is a PLANE CUBIC through the origin.  Hence
\\ Y = {Rpq=0} subset A^2_(p,q) x A^1_c is an elliptic surface with the
\\ zero-section {p=q=0}.  Find its singular fibres and its Weierstrass form.
u='u; s='s; t='t; p='p; q='q; c='c;
A(w0) = (2*w0^5+4*w0^4+6*w0^3+8*w0^2+10*w0+5)/(2*(w0+1)^2);
{
bST = (A(t)-A(s))/(s^2-t^2);
aST = A(s) - bST*(1-s^2);
R = divrem(numerator(A(u) - aST - bST*(1-u^2)), (u-s)*(u-t))[1];
R = R/content(R);
Rpq = subst(subst(R,u,c), t, p-s) % (s^2-p*s+q);
Rpq = subst(Rpq, s, 0); Rpq = Rpq/content(Rpq);

\\ --- put the cubic in Weierstrass form.  Origin (0,0) is on every fibre.
\\ Homogeneous parts:
L  = polcoef(polcoef(Rpq,1,p),0,q)*p + polcoef(polcoef(Rpq,0,p),1,q)*q;
print("linear part  L1 = ", L, "   (factored: ", factor(L), ")");
Q2 = 0; for(i=0,3,for(j=0,3, if(i+j==2, Q2 += polcoef(polcoef(Rpq,i,p),j,q)*p^i*q^j)));
C3 = 0; for(i=0,3,for(j=0,3, if(i+j==3, C3 += polcoef(polcoef(Rpq,i,p),j,q)*p^i*q^j)));
print("quadratic part Q2 = ", Q2);
print("cubic part     C3 = ", C3);
print("");
print("=> origin is a SINGULAR point of the fibre exactly when L1 == 0, i.e. c^2(2c+1)=0:");
print("   c = 0 (double)  and  c = -1/2 .");
print("   c=0 is degenerate for the chart (v=0 <=> root r=1 = the marked point).");
print("   c=-1/2 is the UNIQUE non-degenerate nodal fibre through the section.");
print("");
\\ --- discriminant of the cubic fibration ---
\\ Weierstrass: put p = P, q = Q; use the standard cubic-with-rational-point reduction
\\ by brute force: intersect with lines through origin, p = k*q.
Par = subst(Rpq, p, k*q);
\\ Par = q*(L1coeff) + q^2*(...) + q^3*(...)
c1 = polcoef(Par,1,q); c2 = polcoef(Par,2,q); c3 = polcoef(Par,3,q);
print("line p = k q meets the cubic where  c3 q^2 + c2 q + c1 = 0 with");
print("  c1 = ", c1);
print("  c2 = ", c2);
print("  c3 = ", c3);
D = c2^2 - 4*c3*c1;
D = D/content(D);
print("");
print("residual quadratic discriminant in q:   D(k,c) = ", D);
print("  deg_k = ", poldegree(D,k), "   deg_c = ", poldegree(D,c));
print("  factor over Q(c)[k]: ", factor(D));
}
quit
