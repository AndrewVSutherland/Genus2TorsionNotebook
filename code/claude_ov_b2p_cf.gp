\\ Lane 3 (route B2'): exact CF order of D_inf = inf+ - inf- for monic sextics,
\\ plus the [1,1,2,2] parametrizations of the BLP2009 ansatz.
\\
\\ Charts:
\\   BLP:   y^2 = R(x)^2 - 4 c^2 S(x)^2,  R = x^3 - x^2 + a x + b, S = x^2 + d
\\          f = (R-2cS)(R+2cS) = product of two monic cubics.
\\   B2':   [1,1,2,2] locus = each cubic acquires a rational root r1, r2.
\\          Given (c,d,r1,r2) the two linear equations
\\             (R-2cS)(r1) = 0,  (R+2cS)(r2) = 0
\\          solve for (a,b).
\\   MONIC: y^2 = x (x-1) (x^2+p1 x+p2) (x^2+p3 x+p4)   -- the normalized
\\          [1,1,2,2] chart (two rational Weierstrass pts moved to 0 and 1).
\\
\\ Usage:  gp -q code/claude_ov_b2p_cf.gp
\\ Self-test vectors (must pass): f14 -> 14, f18 -> 18, C4corr -> 11, 19044.h.2 -> 11.

x = 'x;

\\ ---------------------------------------------------------------- CF order
\\ polynomial part of sqrt(f) for f monic of degree 6
sqrtpolypart(f) = {
  my(s = x^3, d);
  for(k=1,3,
    d = f - s^2;
    if(poldegree(d) <= 2, break);
    s = s + (polcoef(d, 6-k)/2) * x^(3-k));
  s;
}

\\ exact order of D_inf; returns 0 if no quasi-period within maxsteps
cford(f, maxsteps=60) = {
  my(s, Pi, Qi, ai, Pn, Qn, tot, r);
  if(poldegree(f) != 6 || pollead(f) != 1, return(-1));
  s = sqrtpolypart(f);
  Pi = 0; Qi = 1; tot = 0;
  for(i=0, maxsteps,
    if(Qi == 0, return(0));
    ai = (Pi + s) \ Qi;
    tot += poldegree(ai);
    Pn = ai*Qi - Pi;
    r = f - Pn^2;
    if(r % Qi != 0, return(0));
    Qn = r / Qi;
    Pi = Pn; Qi = Qn;
    if(i >= 1 && poldegree(Qi) <= 0 && Qi != 0, return(tot)));
  0;
}

\\ ------------------------------------------------- BLP / B2' parametrization
\\ given (c,d,r1,r2) return [a,b] solving (R-2cS)(r1)=0, (R+2cS)(r2)=0
b2p_ab(c,d,r1,r2) = {
  my(a, b);
  if(r1 == r2, return(0));
  a = ( -(r1^3 - r2^3) + (1+2*c)*r1^2 + (2*c-1)*r2^2 + 4*c*d ) / (r1 - r2);
  b = -r1^3 + (1+2*c)*r1^2 - a*r1 + 2*c*d;
  [a,b];
}

b2p_f(c,d,r1,r2) = {
  my(v = b2p_ab(c,d,r1,r2), R, S);
  if(v == 0, return(0));
  R = x^3 - x^2 + v[1]*x + v[2];
  S = x^2 + d;
  R^2 - 4*c^2*S^2;
}

blp_f(a,b,c,d) = (x^3-x^2+a*x+b)^2 - 4*c^2*(x^2+d)^2;

monic_f(p1,p2,p3,p4) = x*(x-1)*(x^2+p1*x+p2)*(x^2+p3*x+p4);

\\ normalize a monic sextic with >= 2 rational roots into the MONIC chart:
\\ send two chosen rational roots to 0 and 1 by x -> u*X + v (affine keeps inf+-)
nrr(g) = my(fa=factor(g), s=0); for(j=1,matsize(fa)[1], if(poldegree(fa[j,1])==1, s+=fa[j,2])); s;
ratroots(g) = my(fa=factor(g), L=[]); for(j=1,matsize(fa)[1], if(poldegree(fa[j,1])==1, L=concat(L,[-polcoef(fa[j,1],0)/polcoef(fa[j,1],1)]))); L;

\\ f monic sextic, roots R1,R2 rational: substitute x = (R2-R1)*X + R1, divide by (R2-R1)^6
to_monic_chart(f, R1, R2) = {
  my(u = R2-R1, g);
  g = subst(f, x, u*x + R1) / u^6;
  g;
}

{
print("=== self-test: validated CF vectors ===");
f14 = (x^2+1)*(x^4+5*x^2+4*x+4);
f18 = (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8);
f28 = x^6+2*x^5-5*x^4-14*x^3-3*x^2+24*x+28;
print("  f14 -> ", cford(f14), "  (expect 14)");
print("  f18 -> ", cford(f18), "  (expect 18)");
print("  f28 -> ", cford(f28), "  (expect 7, the cautionary vector)");
if(cford(f14)!=14 || cford(f18)!=18 || cford(f28)!=7, error("CF SELF-TEST FAILED"));
print("  SELF-TEST OK");

print("");
print("=== anchor 1: BLP C4corr  (a,b,c,d) = (1159/81, -277/243, 40/9, 13/27) ===");
a0=1159/81; b0=-277/243; c0=40/9; d0=13/27;
F = blp_f(a0,b0,c0,d0);
print("  f = ", F);
print("  factored: ", factor(F));
print("  CF order = ", cford(F), "   (expect 11)");
\\ rational roots of the two cubics
A = x^3-x^2+a0*x+b0 - 2*c0*(x^2+d0);
B = x^3-x^2+a0*x+b0 + 2*c0*(x^2+d0);
r1 = ratroots(A); r2 = ratroots(B);
print("  rational roots of R-2cS: ", r1, "   of R+2cS: ", r2);
\\ round-trip through the B2' parametrization
if(#r1>0 && #r2>0,
  ab = b2p_ab(c0,d0,r1[1],r2[1]);
  print("  B2' round trip: (a,b) recovered = ", ab, "   original = ", [a0,b0]);
  if(ab != [a0,b0], error("B2' PARAMETRIZATION ROUND TRIP FAILED"));
  print("  round-trip f matches: ", b2p_f(c0,d0,r1[1],r2[1]) == F);
  print("  CF order from parametrization = ", cford(b2p_f(c0,d0,r1[1],r2[1]))));
\\ integral model quoted in the brief
Fint = (x-9)*(x+21)*(x^2-80*x+439)*(x^2+50*x+109);
print("  brief's integral model CF order = ", cford(Fint), " (scaled model, expect 11)");
print("  scaling check: 9^6*subst(F,x,x/9) == Fint ? ", 9^6*subst(F,x,x/9) == Fint);

print("");
print("=== anchor 2: 19044.h.2, MONIC chart (p1,p2,p3,p4) = (-3,8,4,27) ===");
G = monic_f(-3,8,4,27);
print("  f = ", G);
print("  CF order = ", cford(G), "   (expect 11)");
print("  factor: ", factor(G));

print("");
print("=== anchor 1 in the MONIC chart ===");
rr = ratroots(F);
print("  rational roots of C4corr's f: ", rr);
H = to_monic_chart(F, rr[1], rr[2]);
print("  normalized f = ", H);
print("  CF order = ", cford(H));
print("  as x(x-1)*q1*q2: ", factor(H));

print("");
print("=== C4corr MONIC-chart parameters (p1,p2,p3,p4) ===");
fa = factor(H);
for(j=1,matsize(fa)[1], if(poldegree(fa[j,1])==2, print("  quad: ", fa[j,1])));

print("");
print("=== the other normalization (swap the two roots) ===");
H2 = to_monic_chart(F, rr[2], rr[1]);
print("  CF order = ", cford(H2), "  factor: ", factor(H2));
}
quit
