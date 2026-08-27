Q:=Rationals();
K<r,t>:=RationalFunctionField(Q,2);
Px<x>:=PolynomialRing(K);
pv:=r*(r+t)/2; e:=t^2-2*pv*t/r; lambda:=r/t;
A:=r^2-lambda; B:=2*r*pv-2*lambda*(pv+r*t)+2*r*lambda;
C:=pv^2+2*pv*r^2-r^4-r^3*t-r*pv^2/t-lambda*(r^2+e)+2*lambda*(r*pv+r^2*t-3*pv*t+r*t^2);
qq:=A*x^2+B*x+C; f:=qq*(x^4+qq);
fc:=[Coefficient(f,i):i in [0..6]]; f0:=fc[1];f1:=fc[2];f2:=fc[3];f3:=fc[4];f4:=fc[5];f5:=fc[6];f6:=fc[7];
l:=C; j:=f1/(2*l); n:=(f2-j^2)/(2*l);
S<m,w>:=PolynomialRing(K,2);
kap:=f6-m^2; bN:=2*m*n-f5;
Eq4:=3*kap*(f4-n^2-2*m*j)-bN^2;
Eq3:=27*kap^2*(f3-2*m*l-2*n*j)+bN^3;
I:=ideal<S| Eq4, Eq3, w*kap-1>;   // Rabinowitsch: forces kap<>0
G:=GroebnerBasis(I);
printf "GroebnerBasis size=%o\n", #G;
printf "Ideal is (1) [no genuine soln generically]? %o\n", I eq ideal<S|1>;
if I ne ideal<S|1> then
  printf "GB (genuine-locus relations):\n";
  for g in G do printf "  %o\n", g; end for;
end if;
quit;
