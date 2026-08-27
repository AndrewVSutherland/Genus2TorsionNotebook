\\ verify.gp — Lane C structural verifications on the T5 pencil of M(2,2,2,6).
\\ Variables: v = 2u-1, Q = v^2-v+1 (= q(u)), r = rho', X = Q*r.
\\ Claims to verify (all "==" are equalities of square classes in Q(v,r)):
\\  (F1) V'1 == Q*X*(X-v)
\\  (F2) V'2 == Q*(X-Q)*(X-v)
\\  (F3) V'3 == v*X*(X-Q)*(X-1)*(v*X-Q)
\\  (F4) V'4 == ? (compute; candidates v*Q*(X-v)*(v*X-Q) or v*X*(X-v)*(v*X-Q))
\\  (I)  V'1*V'2*V'3*V'4 == (X-1)*(X-v)
\\  (L)  V'1*V'2 == r*(r-1)   [THE LAW: pass V'1&V'2 => rn(rn-rd) = square]
\\  (S)  sigma: v -> v/(v-1), r -> 1-r: action on classes; on S (t-coords) C1<->C2, D4 fixed
\\  (T)  on S: r = X/Q, X = (t^2-v)/(t^2-1); sigma acts as t -> 1/t
default(parisize, "1G");
Av = [1,1,1,2,2];
BS(s,m,n) = [2*s^2-s*n, 2*s^2+s*m-2*s*n-m*n, 2*s^2+s*m-s*n-m*n, -m*n, 4*s^2-4*s*n-m*n];
\\ square class of a rational function: signed squarefree kernel incl. constant
sqcore(f) = {
  my(fa = factor(f), od = 1, cc);
  for(j=1, matsize(fa)[1], if(fa[j,2] % 2, od *= fa[j,1]));
  cc = simplify(f / factorback(fa));
  core(numerator(cc)*denominator(cc)) * od;
}
samecl(f, g) = issquare(f*g);
{
\\ ---------- build V'_k on the T5 pencil in (v, r) ----------
u = ('v+1)/2;
s = u*('r-1); m = 4*'r*u*(u-1); n = 'r-1;
Bp = BS(s, m, n);
Vk = vector(5);
for(k=1,4,
  my(E = Av[k]*Bp[5] - Av[5]*Bp[k], P = 1);
  for(l=1,4, if(l!=k, P *= Bp[l]));
  Vk[k] = E*P);
Q = 'v^2-'v+1;  X = Q*'r;
F = [Q*X*(X-'v), Q*(X-Q)*(X-'v), 'v*X*(X-Q)*(X-1)*('v*X-Q), 0];
print("=== factored normal forms in (v, X=Q*r) ===");
for(k=1,3, printf("V'%d == claimed: %d\n", k, samecl(Vk[k], F[k])));
\\ V'4: compute its class and try candidates
c4a = 'v*Q*(X-'v)*('v*X-Q);
c4b = 'v*X*(X-'v)*('v*X-Q);
printf("V'4 == v*Q*(X-v)*(vX-Q): %d\n", samecl(Vk[4], c4a));
printf("V'4 == v*X*(X-v)*(vX-Q): %d\n", samecl(Vk[4], c4b));
printf("V'4 core = %s\n", Str(sqcore(Vk[4])));
\\ product identity
printf("prod V' == (X-1)(X-v): %d\n", samecl(Vk[1]*Vk[2]*Vk[3]*Vk[4], (X-1)*(X-'v)));
\\ THE LAW
printf("LAW  V'1*V'2 == r*(r-1): %d\n", samecl(Vk[1]*Vk[2], 'r*('r-1)));
printf("     V'3*V'4 core = %s\n", Str(sqcore(Vk[3]*Vk[4])));
printf("     V'1*V'3 core = %s\n", Str(sqcore(Vk[1]*Vk[3])));
printf("     V'2*V'3 core = %s\n", Str(sqcore(Vk[2]*Vk[3])));
printf("     V'1*V'4 core = %s\n", Str(sqcore(Vk[1]*Vk[4])));
printf("     V'2*V'4 core = %s\n", Str(sqcore(Vk[2]*Vk[4])));
\\ ---------- sigma off S: class transport ----------
print("=== sigma: v -> v/(v-1), r -> 1-r (X -> (Q-X)/(v-1)^2) ===");
sig(f) = substvec(f, ['v,'r], ['v/('v-1), 1-'r]);
for(k=1,4, printf("class(V'%d o sigma) core = %s\n", k, Str(sqcore(sig(Vk[k])))));
for(k=1,4, for(l=1,4,
  if(samecl(sig(Vk[k]), Vk[l]), printf("  V'%d o sigma == V'%d (globally)\n", k, l))));
\\ ---------- on S: t-parameterization ----------
print("=== S-surface: X = (t^2-v)/(t^2-1), i.e. r = X/Q ===");
Xt = ('t^2-'v)/('t^2-1);
rt = Xt/Q;
VS = vector(4, k, substvec(Vk[k], ['r], [rt]));
C1 = Q*(1-'v)*('t^2-'v);
C2 = Q*('v*'t^2-'v+1);
D4 = 'v*(('v-1)*'t^2+1);
printf("on S: V'1 == C1 = Q(1-v)(t^2-v): %d\n", samecl(VS[1], C1));
printf("on S: V'2 == C2 = Q(v t^2-v+1): %d\n", samecl(VS[2], C2));
printf("on S: V'4 == D4 = v((v-1)t^2+1): %d\n", samecl(VS[4], D4));
printf("on S: V'3 == C1*C2*D4: %d\n", samecl(VS[3], C1*C2*D4));
printf("on S: V'3 core = %s\n", Str(sqcore(VS[3])));
\\ sigma in (v,t): v -> v/(v-1), t -> 1/t
sigS(f) = substvec(f, ['v,'t], ['v/('v-1), 1/'t]);
printf("sigma(C1) == C2: %d ; sigma(C2) == C1: %d ; sigma(D4) == D4: %d\n",
  samecl(sigS(C1), C2), samecl(sigS(C2), C1), samecl(sigS(D4), D4));
\\ consistency: sigma(X(t,v)) equals X(1/t, v/(v-1))?
X2a = substvec(Xt, ['v,'t], ['v/('v-1), 1/'t]);
X2b = (Q - Xt)/('v-1)^2;   \\ sigma on X derived from r->1-r
printf("t->1/t consistency on S: %d\n", X2a == X2b);
\\ hit coordinates in (v,t):
vh = -121/24; Xh = subst(subst(Xt,'v,vh),'t,437/1013);
printf("hit check: X(v=-121/24, t=437/1013) = %s (expect -887929/138240): %d\n",
  Str(Xh), Xh == -887929/138240);
\\ values of C1,C2,D4 at the hit (should all be squares):
for(i=1,3, my(Ci=[C1,C2,D4][i], val=subst(subst(Ci,'v,vh),'t,437/1013));
  printf("  %s at hit = %s, issquare=%d\n", ["C1","C2","D4"][i], Str(val), issquare(val)));
}
quit
