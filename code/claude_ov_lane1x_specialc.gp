\\ claude_ov_lane1x_specialc.gp -- Lane 1: EXACTLY which u = c slices of the
\\ contact-7 three-root surface fail to have genus 4.
\\
\\ C_c = {R(s,t,c)=0} in A^2 is a (3,3) curve, arithmetic genus 4.  Under
\\ (p,q)=(s+t,st) it is the DOUBLE COVER of the plane cubic E_c = {Rpq(p,q,c)=0}
\\ branched at the 6 points of E_c cap {p^2=4q}.  Hence genus(C_c) = 4 unless
\\   (i)  E_c is singular            [disc of the cubic vanishes], or
\\   (ii) the branch divisor is non-reduced  [E_c tangent to the conic p^2=4q].
\\ Both are polynomial conditions in c; compute them and factor over Q.
u='u; s='s; t='t; p='p; q='q; c='c; k='k;
A(w0) = (2*w0^5+4*w0^4+6*w0^3+8*w0^2+10*w0+5)/(2*(w0+1)^2);
ratroots(F) = my(fa=factor(F), r=List()); for(i=1,#fa[,1], if(poldegree(fa[i,1],c)==1, listput(r, -polcoef(fa[i,1],0,c)/polcoef(fa[i,1],1,c)))); Vec(r);
{
bST = (A(t)-A(s))/(s^2-t^2);
aST = A(s) - bST*(1-s^2);
R = divrem(numerator(A(u) - aST - bST*(1-u^2)), (u-s)*(u-t))[1];
R = R/content(R);
Rpq = subst(subst(R,u,c), t, p-s) % (s^2-p*s+q);
Rpq = subst(Rpq, s, 0); Rpq = Rpq/content(Rpq);
print("Rpq(p,q,c) = ", Rpq);
print("");

\\ ---------- (i) singular fibres of the cubic pencil ----------
\\ E_c passes through (0,0) for every c; project from the origin: p = k q.
Par = subst(Rpq, p, k*q)/q;      \\ = c1 + c2 q + c3 q^2
cc1 = polcoef(Par,0,q); cc2 = polcoef(Par,1,q); cc3 = polcoef(Par,2,q);
Dk  = cc2^2 - 4*cc3*cc1;         \\ E_c is y^2 = Dk(k), a cubic in k
print("E_c as a double cover of the k-line:  y^2 = Dk(k,c), deg_k = ", poldegree(Dk,k));
dis = poldisc(Dk, k);
dis = dis/content(dis);
print("");
print("disc_k(Dk) as a polynomial in c:  degree ", poldegree(dis,c));
fa = factor(dis);
print("  factorisation:");
for(i=1,#fa[,1], print("    (", fa[i,1], ")^", fa[i,2]));
print("  RATIONAL roots c (singular fibres of the cubic pencil): ", ratroots(dis));

\\ ---------- (ii) tangency of E_c with the conic p^2 = 4q ----------
\\ substitute q = p^2/4 into Rpq; the 6 branch points are the roots in p
Br = numerator(subst(Rpq, q, p^2/4));
Br = Br/content(Br);
print("");
print("branch polynomial (E_c cap {p^2=4q}) in p:  degree ", poldegree(Br,p));
print("  Br = ", Br);
dbr = poldisc(Br, p);
dbr = dbr/content(dbr);
print("");
print("disc_p(Br) as a polynomial in c:  degree ", poldegree(dbr,c));
fb = factor(dbr);
print("  factorisation:");
for(i=1,#fb[,1], print("    (", fb[i,1], ")^", fb[i,2]));
print("  RATIONAL roots c (tangency / non-reduced branch divisor): ", ratroots(dbr));
print("");
print("  leading coefficient of Br in p (drop of degree at infinity): ",
      factor(polcoef(Br, poldegree(Br,p), p)));
}
quit
