\\ v2qf.gp — DS-v2 quintic-cofactor quadratic-factor incidence via RESULTANTS
\\ (replaces the Magma 3-variable Groebner formulation that hit the memory cap).
\\ g_u has the rational root -u/2; q5 = g_u/(x+u/2).  Incidence:
\\ x^2+A*x+B | q5  <=>  e1(u,A,B) = e0(u,A,B) = 0 (division remainder).
\\ Eliminate B: RB = Res_B(e1,e0) in Q[u,A]; factor; report component shapes.
default(parisizemax, 8*10^9);
x='x; u='u; A='A; B='B;
{
den = u^5+8;
g = x^6 - (16*u/den)*x^5 \
 - ((u^15+40*u^10+128*u^5+512)/(2*u^3*den^2))*x^4 \
 + ((6*u^15+224*u^10+1536*u^5+3072)/(u^2*den^3))*x^3 \
 + ((u^30+112*u^25+2880*u^20+25600*u^15+106496*u^10+262144*u^5+262144)/(16*u^6*den^4))*x^2 \
 - ((u^25+80*u^20+1664*u^15+12288*u^10+36864*u^5+32768)/(2*u^5*den^4))*x \
 - (u^25+46*u^20+736*u^15+5248*u^10+18432*u^5+24576)/(2*u^4*den^4);
dv = divrem(g, x + u/2);
if(dv[2] != 0, print("ERROR: -u/2 not a root"); quit);
q5 = dv[1];
dd = lcm(vector(5, j, denominator(polcoef(q5, j-1, x))));
q5i = dd*q5;
print("q5 cleared; coefficient u-degrees: ", vector(6, j, poldegree(polcoef(q5i, j-1, x), u)));
dv2 = divrem(q5i, x^2 + A*x + B);
r = dv2[2];
e1 = polcoef(r, 1, x); e0 = polcoef(r, 0, x);
print("e1 degs (u,A,B): ", [poldegree(e1,u), poldegree(e1,A), poldegree(e1,B)]);
print("e0 degs (u,A,B): ", [poldegree(e0,u), poldegree(e0,A), poldegree(e0,B)]);
RB = polresultant(e1, e0, B);
RB = RB/content(RB);
print("Res_B computed: degs (u,A) = ", [poldegree(RB,u), poldegree(RB,A)]);
fw = factor(RB);
print("factors (mult, deg_u, deg_A):");
for(k=1, matsize(fw)[1],
  gk = fw[k,1];
  if(poldegree(gk,u) == 0 && poldegree(gk,A) == 0, next);
  print("  [", fw[k,2], "] deg_u=", poldegree(gk,u), " deg_A=", poldegree(gk,A));
  if(poldegree(gk,A) <= 2 && poldegree(gk,u) <= 40, print("    ", gk)));
print("V2QF_RESULTANT_DONE");
}
quit
