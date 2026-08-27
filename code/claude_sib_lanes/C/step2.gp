\\ step2.gp — sigma action + S-surface forms, using the VERIFIED factored classes:
\\ V'1 == Q*X*(X-v), V'2 == Q*(X-Q)*(X-v), V'3 == v*Q*X*(X-1)*(X-Q)*(vX-Q),
\\ V'4 == v*Q*(X-v)*(vX-Q)   (all mod squares; X = Q*r, Q = v^2-v+1)
default(parisize, "512M");
samecl(f, g) = issquare(f*g);
sig(f) = substvec(f, ['v,'r], ['v/('v-1), 1-'r]);
sigS(f) = substvec(f, ['v,'t], ['v/('v-1), 1/'t]);
sqcore(f) = {
  my(fa = factor(f), od = 1, cc);
  for(j=1, matsize(fa)[1], if(fa[j,2] % 2, od *= fa[j,1]));
  cc = simplify(f / factorback(fa));
  core(numerator(cc)*denominator(cc)) * od;
}
{
Q = 'v^2-'v+1;  X = Q*'r;
F = [Q*X*(X-'v), Q*(X-Q)*(X-'v), 'v*Q*X*(X-1)*(X-Q)*('v*X-Q), 'v*Q*(X-'v)*('v*X-Q)];
\\ verify V'3 corrected form against product identity: F1*F2*F3*F4 == (X-1)(X-v)
printf("F-product == (X-1)(X-v): %d\n", samecl(F[1]*F[2]*F[3]*F[4], (X-1)*(X-'v)));
\\ sigma: v -> v/(v-1), r -> 1-r
print("=== sigma class transport (global, off S) ===");
for(k=1,4, my(g = sig(F[k]), hit=0);
  for(l=1,4, if(samecl(g, F[l]), printf("V'%d o sigma == V'%d\n", k, l); hit=1));
  if(!hit, printf("V'%d o sigma == (new class) core = %s\n", k, Str(sqcore(g)))));
\\ pairwise products transported
print("=== sigma on products ===");
PR = [F[1]*F[2], F[3]*F[4], F[1]*F[3], F[1]*F[4], F[2]*F[3], F[2]*F[4]];
PRn = ["12","34","13","14","23","24"];
for(i=1,6, my(g = sig(PR[i]), hit=0);
  for(j=1,6, if(samecl(g, PR[j]), printf("V'%s o sigma == V'%s\n", PRn[i], PRn[j]); hit=1));
  if(!hit, printf("V'%s o sigma: new class\n", PRn[i])));
\\ ---------- on S ----------
print("=== S-surface (t-coords) ===");
Xt = ('t^2-'v)/('t^2-1);
FS = vector(4, k, substvec(F[k], ['r], [Xt/Q]));
C1 = Q*(1-'v)*('t^2-'v);
C2 = Q*('v*'t^2-'v+1);
D4 = 'v*(('v-1)*'t^2+1);
printf("on S: V'1 == C1: %d ; V'2 == C2: %d ; V'4 == D4: %d ; V'3 == C1*C2*D4: %d\n",
  samecl(FS[1], C1), samecl(FS[2], C2), samecl(FS[4], D4), samecl(FS[3], C1*C2*D4));
printf("sigma(C1)==C2: %d ; sigma(C2)==C1: %d ; sigma(D4)==D4: %d\n",
  samecl(sigS(C1), C2), samecl(sigS(C2), C1), samecl(sigS(D4), D4));
\\ t->1/t consistency with r->1-r on S
X2a = substvec(Xt, ['v,'t], ['v/('v-1), 1/'t]);
X2b = substvec((Q - X)/('v-1)^2, ['r], [Xt/Q]);
printf("sigma on S is t->1/t: %d\n", X2a - X2b == 0);
\\ hit coordinates
vh = -121/24; th = 437/1013;
Xh = subst(subst(Xt,'v,vh),'t,th);
printf("hit: X(v=-121/24,t=437/1013) = -887929/138240: %d\n", Xh == -887929/138240);
for(i=1,3, my(Ci=[C1,C2,D4][i], val=subst(subst(Ci,'v,vh),'t,th));
  printf("  %s at hit: issquare=%d\n", ["C1","C2","D4"][i], issquare(val)));
\\ ---------- laws recap on S ----------
\\ r(r-1) on S:
rt = Xt/Q;
printf("on S: r*(r-1) == -C1*C2 (i.e. == v-1 ... ): core(r(r-1)*C1*C2) = %s\n",
  Str(sqcore(rt*(rt-1)*C1*C2)));
\\ D4 as conic in t with point t=1: parameterization t = ((k-v)^2-v)/(k^2-v^2+v)
tk = (('k-'v)^2-'v)/('k^2-'v^2+'v);
D4k = substvec(D4, ['t], [tk]);
printf("D4 on k-param is a square: %d\n", issquare(D4k));
C1k = substvec(C1, ['t], [tk]);
C2k = substvec(C2, ['t], [tk]);
print("C1 pulled to (v,k): core = ", Str(sqcore(C1k)));
print("C2 pulled to (v,k): core = ", Str(sqcore(C2k)));
print("C1*C2 pulled: core = ", Str(sqcore(C1k*C2k)));
print("factor C1k core: ", factor(sqcore(C1k)));
print("factor C2k core: ", factor(sqcore(C2k)));
}
quit
