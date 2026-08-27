\\ claude_c7_relations.gp -- the three universal relations governing the rational
\\ Weierstrass points of the contact-7 chart, and their verification.
\\
\\ Chart (notes/contact7_family.md):  h = 1 - (7/2)x + a x^2 + b x^3,
\\                                    f = (h^2 + (x-1)^7)/x^2  (monic quintic),
\\ marked class [P - inf] of order 7 at P = (1, h(1)).
\\ A rational root r of f forces h(r)^2 = (1-r)^7, hence 1-r = v^2 with v rational
\\ and h(r) = v^7.  The five v_i are the roots of
\\    Q(v) = (v+1)^2 (c4 v^2 + c0) + v^5 - v^3 - v^2/2,   c4 = b+2, c0 = 5/2-(a+b).
\\ Only two parameters enter, so the elementary symmetric functions e1..e5 of the
\\ five roots satisfy three universal relations, verified below for random (a,b):
\\    (R1) e2 + 2 e1 + 1     = 0
\\    (R2) e4 + 2 e5         = 0        (i.e. sum 1/v_i = -2)
\\    (R3) 2 e3 - 2 e1 - 2 e5 - 1 = 0
\\ k rational v_i  <=>  quintic factor type with k rational roots  <=>  2-rank k-1
\\ (k<5) or 4 (k=5); with the order-7 class:
\\    k=1 [14] | k=2 [2,14] | k=3 [2,2,14] (realized 2026-07-23) | k=5 [2,2,2,14], order 112.
G(u) = -(u^5 - u^3 - u^2/2)/(u+1)^2;
Qpoly(c4,c0) = v^5 + c4*v^4 + (2*c4-1)*v^3 + (c4+c0-1/2)*v^2 + 2*c0*v + c0;
getc(v1,v2) = my(c4,c0); c4=(G(v1)-G(v2))/(v1^2-v2^2); c0=G(v1)-c4*v1^2; [c4,c0];
\\ --- 1. relations hold identically ---
for(t=1,5, my(aa,bb,c4,c0,Q,e1,e2,e3,e4,e5); aa=random(1000)/(1+random(97)); bb=random(1000)/(1+random(89)); c4=bb+2; c0=5/2-(aa+bb); Q=Qpoly(c4,c0); e1=-polcoef(Q,4); e2=polcoef(Q,3); e3=-polcoef(Q,2); e4=polcoef(Q,1); e5=-polcoef(Q,0); print("(a,b)=(",aa,",",bb,")  R1=",e2+2*e1+1,"  R2=",e4+2*e5,"  R3=",2*e3-2*e1-2*e5-1));
\\ --- 2. the six originally recorded three-root hits, plus the five found by the harvester ---
hits = [[-3,-3/4,-3/5],[-10,-1/2,-10/7],[-5,-15/8,-15/22],[-1/2,-15/8,-15/19],[-4/9,4/17,-4/25],[4/17,-5/18,-10/49],[-511/61,-511/625,-1/2],[-165/41,-33/16,-165/289],[-164/297,-1/2,164/361],[-17/50,-34/189,34/121],[-1/2,-13/49,13/50]];
for(i=1,length(hits), my(h,c,Q,fa,degs); h=hits[i]; c=getc(h[1],h[2]); Q=Qpoly(c[1],c[2]); fa=factor(Q); degs=vector(length(fa[,1]),j,poldegree(fa[j,1])); print("hit ",i," v=",h," (c4,c0)=",c," v-quintic factor degrees=",degs));
\\ --- 3. local test for the SPLITALL (order-112) locus: how many (c4,c0) mod p split completely ---
forprime(p=5,60, if(p==7,,my(cnt=0,i2=Mod(1,p)/2); for(A=0,p-1, for(B=0,p-1, my(Q,fa,ok); Q=Mod(1,p)*(v^5+A*v^4+(2*A-1)*v^3+(A+B-1/2)*v^2+2*B*v+B); if(poldisc(Q)!=Mod(0,p), fa=factormod(Q,p); ok=1; for(j=1,length(fa[,1]), if(poldegree(fa[j,1])>1,ok=0)); if(ok,cnt++)))); print("p=",p," split-all (c4,c0) pairs: ",cnt," / ",p^2," (naive heuristic p^2/120 = ",p^2/120.,")")));
quit
