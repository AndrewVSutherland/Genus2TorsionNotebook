SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned log_file then
 log_file:="results/target_22224_W1_exact_survivor.log";
end if;
SetLogFile(log_file:Overwrite:=true);
fixed:=[Q!-2254,-2162,2303]; d0:=Q!4900;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[3]+d0*T^2)/(fixed[3]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
W:=E![Q!-1,Q!0,Q!1];
G1:=E![Q!577,Q!-14552,Q!1];
G2:=E![Q!17249/49,Q!-2560104/343,Q!1];
G3:=E![Q!1873/49,Q!407464/343,Q!1];
ep:=2*G1-G2+G3+W;
cp:=Einv(minmapinv(ep)); assert cp[3] ne 0;
tt:=Q!(cp[1]/cp[3]); dd:=d0*tt^2;
flags:=[]; for z in fixed do ok,s:=IsSquare((z+dd)/(z+d0)); Append(~flags,ok); end for;
print "W1_EXACT_SURVIVOR","coeffs",<2,-1,1,2>,"t",tt,"square_flags",flags;
if &and flags then
 den:=LCM([Denominator(z):z in fixed cat [dd]]);
 v:=[Z!(den*z):z in fixed cat [dd]]; g:=GCD([Abs(z):z in v]); v:=[z div g:z in v];
 P<x>:=PolynomialRing(Q); f:=x*&*[x+(Q!z)^2:z in v];
 JG,jm:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
 print "W1_EXACT_TORSION",v,Invariants(JG),#JG,f;
end if;
UnsetLogFile(); quit;
