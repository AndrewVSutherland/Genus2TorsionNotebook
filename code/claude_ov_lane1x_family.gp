\\ claude_ov_lane1_family.gp -- Lane 1: the INFINITE [2,2,14] family on the
\\ u = -1/2 slice of the contact-7 three-root surface.
\\
\\   R(s,t,-1/2)=0  <=>  2q^3+3pq^2+2pq-p^3 = 0   (p=s+t, q=st): NODAL cubic
\\   p = 2m^2/((m+1)^2(m-2)),  q = 2m/((m+1)^2(m-2))
\\   s,t rational  <=>  w^2 = -m^4+6m^2+4m = -m(m+2)(m^2-2m-2)
\\   s = (m^2+w)/((m+1)^2(m-2)),  t = (m^2-w)/((m+1)^2(m-2))
\\   E: Y^2 = X^3+6X^2-16,  X = 4/m, Y = 4w/m^2.   E(Q) = Z<(-4,4)> + Z/2<(-2,0)>
\\ Output: N family members as (s,t,-1/2) with (a,b) and the quintic f.
if(type(NMAX) != "t_INT", NMAX = 14);
E = ellinit([0,6,0,0,-16]);
P = [-4,4]; T = [-2,0];
Afun(w0) = (2*w0^5+4*w0^4+6*w0^3+8*w0^2+10*w0+5)/(2*(w0+1)^2);
x='x;
{
print("E: ", E.disc, "  conductor ", ellglobalred(E)[1], "  rank 1 gen P=",P,"  T=",T);
print("");
out = List();
for(n=-NMAX, NMAX,
 for(tw=0,1,
  my(Q0, m, w, den, s, t, b, a, h, F7, dq, f, fa, degs, key);
  Q0 = ellmul(E,P,n);
  if(tw==1, Q0 = elladd(E,Q0,T));
  if(#Q0 < 2, next);                       \\ point at infinity
  if(Q0[1] == 0, next);
  m = 4/Q0[1];
  if(m == 0 || m == -1 || m == 2, next);
  w = Q0[2]*m^2/4;
  if(w^2 != -m^4+6*m^2+4*m, print("W-FAIL n=",n," tw=",tw); next);
  den = (m+1)^2*(m-2);
  s = (m^2+w)/den; t = (m^2-w)/den;
  \\ admissibility: the three v's must be distinct, none 0,-1, and no v_i=-v_j
  vs = [s,t,-1/2];
  bad = 0;
  for(i=1,3, if(vs[i]==0 || vs[i]==-1, bad=1));
  for(i=1,3, for(j=i+1,3, if(vs[i]==vs[j] || vs[i]==-vs[j], bad=1)));
  if(bad, next);
  b = (Afun(t)-Afun(s))/(s^2-t^2);
  a = Afun(s) - b*(1-s^2);
  h = 1 - 7/2*x + a*x^2 + b*x^3;
  F7 = h^2 + (x-1)^7;
  dq = divrem(F7, x^2);
  if(dq[2] != 0, print("XDIV-FAIL n=",n); next);
  f = dq[1];
  if(poldegree(f) != 5 || poldisc(f) == 0, next);
  fa = factor(f);
  degs = vecsort(vector(#fa[,1], k, poldegree(fa[k,1])));
  \\ dedupe on the unordered v-set
  key = vecsort([s,t,-1/2]);
  listput(out, [n, tw, m, s, t, a, b, degs, f]);
 ));
\\ dedupe
seen = List(); res = List();
for(k=1,#out, my(key = vecsort([out[k][4],out[k][5]])); my(dup=0);
  for(j=1,#seen, if(seen[j]==key, dup=1));
  if(!dup, listput(seen,key); listput(res,out[k])));
print("distinct family members generated: ", #res);
print("");
nOK = 0;
for(k=1,#res, my(r=res[k]);
  print("n=",r[1]," tw=",r[2]," m=",r[3]);
  print("   (s,t,u) = (",r[4],", ",r[5],", -1/2)");
  print("   (a,b)   = (",r[6],", ",r[7],")");
  print("   f-factor-degrees = ",r[8], if(r[8]==[1,1,1,2]," <-- [1,1,1,2]  2-rank 3"," <-- OTHER"));
  if(r[8]==[1,1,1,2], nOK++);
  print("   f = ", r[9]);
);
print("");
print("members with factor type [1,1,1,2]: ", nOK, " / ", #res);
\\ machine-readable dump for the Magma verifier
write("data/claude_ov_lane1_family.txt", "");
for(k=1,#res, my(r=res[k]);
  write("data/claude_ov_lane1_family.txt", r[1],",",r[2],",",r[3],",",r[4],",",r[5],",",r[6],",",r[7]));
print("wrote data/claude_ov_lane1_family.txt");
}
quit
