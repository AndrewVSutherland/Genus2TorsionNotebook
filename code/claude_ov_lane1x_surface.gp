\\ claude_ov_lane1_surface.gp -- Lane 1 (overnight 2026-07-25):
\\ rebuild the contact-7 three-root surface R(s,t,u), re-verify the recorded
\\ form and all ELEVEN known rational points, and probe the u = -1/2 slice.
u='u; s='s; t='t;
A(w) = (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2);
{
bST = (A(t)-A(s))/(s^2-t^2);
aST = A(s) - bST*(1-s^2);
N = numerator(A(u) - aST - bST*(1-u^2));
dv = divrem(N, (u-s)*(u-t));
if(dv[2] != 0, print("UNEXPECTED remainder"); quit);
R = dv[1];  R = R/content(R);
print("== R rebuilt ==");
print("tridegree = (", poldegree(R,s),",",poldegree(R,t),",",poldegree(R,u),")");
print("R = ", R);
print("S3 checks: cyc=", substvec(R,[s,t,u],[t,u,s])==R, " transp=", substvec(R,[s,t,u],[t,s,u])==R);

\\ --- recorded e-form check ---
\\ R_rec = 2(e2-2)e3^2 + (4e2^2+4e1e2-2e1-4e2-1)e3 + e2(2(e1+e2)^2+e1)
e1 = s+t+u; e2 = s*t+s*u+t*u; e3 = s*t*u;
Rrec = 2*(e2-2)*e3^2 + (4*e2^2+4*e1*e2-2*e1-4*e2-1)*e3 + e2*(2*(e1+e2)^2+e1);
print("recorded-e-form matches R:  ", Rrec == R, "   (ratio ", if(Rrec!=0, R/Rrec, "n/a"), ")");
\\ second recorded form
Rrec2 = 2*e2*(e1+e2+e3)^2 + (s+t)*(t+u)*(u+s) - 2*e3*(e1+2*e2+2*e3);
print("recorded-e-form-2 matches R:", Rrec2 == R);

\\ --- eleven known rational points ---
hits = [[-10,-10/7,-1/2],[-5,-15/8,-15/22],[-3,-3/4,-3/5],[-15/8,-15/19,-1/2],
        [-5/18,-10/49,4/17],[-4/9,-4/25,4/17],[-511/61,-511/625,-1/2],
        [-165/41,-33/16,-165/289],[-164/297,-1/2,164/361],[-17/50,-34/189,34/121],
        [-1/2,-13/49,13/50]];
print("== eleven points ==");
for(k=1,#hits, print("  P",k," = ",hits[k],"   R = ", substvec(R,[s,t,u],hits[k])));

\\ --- the u = -1/2 slice ---
S = subst(R, u, -1/2);  S = S/content(S);
print("== slice u = -1/2 ==");
print("bidegree (s,t) = (", poldegree(S,s),",",poldegree(S,t),")");
print("S = ", S);
print("factor(S) = ", factor(S));
\\ points on the slice
sl = [[-10,-10/7],[-15/8,-15/19],[-511/61,-511/625],[-164/297,164/361],[-13/49,13/50]];
for(k=1,#sl, print("  slice pt ",sl[k]," -> ", substvec(S,[s,t],sl[k])));

\\ --- other fixed-coordinate slices seen among the 11 points ---
print("== other repeated coordinate values ==");
vals = [];
for(k=1,#hits, for(j=1,3, vals = concat(vals,[hits[k][j]])));
vals = vecsort(vals);
cnt = Map();
for(k=1,#vals, my(c); c = if(mapisdefined(cnt,vals[k]), mapget(cnt,vals[k]), 0); mapput(cnt,vals[k],c+1));
for(k=1,#vals, if(k==1 || vals[k]!=vals[k-1], my(c=mapget(cnt,vals[k])); if(c>1, print("  value ",vals[k]," occurs ",c," times"))));
}
quit
