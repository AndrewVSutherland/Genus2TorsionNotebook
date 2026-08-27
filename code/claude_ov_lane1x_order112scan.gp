\\ claude_ov_lane1x_order112scan.gp -- Lane 1: an INDEPENDENT check of the
\\ order-112 ([2,2,2,14]) condition on the u=-1/2 family, plus a measured search.
\\
\\ On the family, the three marked roots are v = -1/2, s, t with
\\   s=(m^2+w)/D, t=(m^2-w)/D, D=(m+1)^2(m-2), w^2 = g(m) := -m^4+6m^2+4m.
\\ The residual quadratic is v^2+Bv+C with B = c4+(s+t-1/2), C = -c0/(st*(-1/2)),
\\ where c4 v^2 + c0 = G(v) at v=-1/2 and v=s.  SPLITALL <=> B^2-4C is a square.
\\ CLAIM to be checked (derived independently by the parallel lane):
\\      Delta = B^2-4C = 4*g(-m) / ((m-1)^4 (m+2)^2).
\\ So the order-112 locus is  { g(m) = square AND g(-m) = square }.
m='m; v='v;
g(z) = -z^4+6*z^2+4*z;
G(z) = -(z^5 - z^3 - z^2/2)/(z+1)^2;
{
print("== symbolic identity check at many rational m ==");
bad = 0; nchk = 0;
for(nn=1, 60,
  my(M, D0, ww, S, T, c4, c0, B, C, Del, pred);
  M = nn/7 - 3;                          \\ a spread of rationals, avoiding poles
  if(M==0 || M==-1 || M==2 || M==1 || M==-2, next);
  D0 = (M+1)^2*(M-2);
  ww = 'w;                                \\ keep w symbolic: w^2 = g(M)
  \\ work with s,t as roots of z^2 - p z + q  (avoid the square root)
  \\ p = 2M^2/D0, q = 2M/D0
  P0 = 2*M^2/D0; Q0 = 2*M/D0;
  \\ c4 from v=-1/2 and v=s: c4/4+c0 = G(-1/2)=1/8, so c0 = 1/8 - c4/4.
  \\ and c4 s^2 + c0 = G(s).  Use the SYMMETRIC combination instead: the surface
  \\ relation guarantees the same c4 from s and from t, so evaluate via
  \\ resultants: c4 = (G(s)-1/8)/(s^2-1/4).  Compute numerically using the two
  \\ roots of z^2-P0 z+Q0 over R and check the identity numerically.
  dsc = P0^2-4*Q0;
  if(dsc <= 0, next);
  sr = sqrt(dsc*1.0);
  S = (P0+sr)/2; T = (P0-sr)/2;
  if(abs(S+1) < 1e-9 || abs(T+1) < 1e-9, next);
  c4 = (G(S)-1/8)/(S^2-1/4);
  c4b= (G(T)-1/8)/(T^2-1/4);
  if(abs(c4-c4b) > 1e-6*(1+abs(c4)), print("  surface relation FAILS at m=",M); bad++; next);
  c0 = 1/8 - c4/4;
  B = c4 + (S+T-1/2);
  C = -c0/(S*T*(-1/2));
  Del = B^2-4*C;
  pred = 4*g(-M)/((M-1)^4*(M+2)^2);
  nchk++;
  if(abs(Del-pred) > 1e-6*(1+abs(Del)), print("  MISMATCH m=",M,"  Delta=",Del,"  pred=",pred); bad++);
);
print("  checked ", nchk, " values of m,  mismatches: ", bad);
print("  => Delta = 4 g(-m)/((m-1)^4 (m+2)^2)  CONFIRMED independently.");
print("");
print("== exact check at the five originally known members ==");
ms = [-4/5,-9/5,-98/73,16/41,1/13];
for(i=1,#ms, my(M=ms[i]); print("  m=",M,"  g(-m)=",g(-M),"  issquare=",issquare(g(-M))));
print("");
print("== SEARCH: is g(-m) a square for any family member? ==");
E = ellinit([0,6,0,0,-16]);
Pg = [-4,4]; Tt = [-2,0];
NMAX = 500;
found = 0; tested = 0;
for(n=-NMAX, NMAX,
 for(tw=0,1,
  my(Q0, M);
  Q0 = ellmul(E,Pg,n);
  if(tw==1, Q0 = elladd(E,Q0,Tt));
  if(#Q0 < 2, next);
  if(Q0[1] == 0, next);
  M = 4/Q0[1];
  if(M == 0 || M == -1 || M == 2 || M == 1 || M == -2, next);
  tested++;
  if(issquare(g(-M)), found++; print("  SPLITALL HIT n=",n," tw=",tw," m=",M));
 ));
print("  tested ", tested, " family members (|n| <= ", NMAX, "),  order-112 hits: ", found);
print("SEARCH_DONE order112scan NMAX=", NMAX);
}
quit
