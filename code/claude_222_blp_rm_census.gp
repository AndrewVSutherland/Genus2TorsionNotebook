\\ RM-signature census on BLP2009 Table 1: real-Weil-subfield disc kernel D_p
\\ across good primes; constant D_p => RM(sqrt D) signature, varying => End=Z signature.
x = 'x;
rows = [[-101/48,-61/48,1/4,-5/12], [473/147,-4013/343,6/7,207/49], [8/49,-134/49,3/7,47/49], \
 [1159/81,-277/243,40/9,13/27], [-1/13,-191/2197,8/13,15/169], [-28/169,103/2197,3/13,-4/169], \
 [594/1805,13348/34295,8/19,-64/361], [208/867,1338/4913,5/17,-39/289], [415/1089,-2207/1089,8/33,119/121], \
 [4989/2500,-13599/12500,27/50,-81/250], [-3,59,4,-7], [-163/1215,-367/3645,2/3,13/243], \
 [-13/18,71/6,5/3,-13/3], [-2287/27,-1171/9,10/3,-323/3], [121/147,-141/343,2/7,15/49], \
 [-1494/847,19480/9317,2/11,-256/121], [125/121,-223/1331,6/11,29/121], [187/361,-649/6859,6/19,23/361]];
names = ["C1","C2","C3","C4corr","C5=X0(23)","C6","C7","C8","C9","C10","Ct1","Ct2","Ct3","Ct4","Ct7","Ct8","Ct9","XX"];
\\ note: row list = C1..C10 with C4 CORRECTED, then Ct1..Ct4, Ct5, Ct7, Ct8, Ct9 (Ct6=Ct3 dropped)
names = ["C1","C2","C3","C4corr","C5=X023","C6","C7","C8","C9","C10","Ct1","Ct2","Ct3","Ct4","Ct5","Ct7","Ct8","Ct9"];
{
\\ control: 1192.a (End=Z [22] curve)
ctrl = 4*(x^3-2*x^2-x+1) + (x^3+x)^2;
curves = concat([[ctrl]], vector(#rows, i, [ (x^3-x^2+rows[i][1]*x+rows[i][2])^2 - 4*rows[i][3]^2*(x^2+rows[i][4])^2 ]));
cnames = concat(["1192a-CTRL"], names);
for(i=1, #curves,
 f = curves[i][1];
 Ds = Set();
 forprime(p=11, 97,
  chi = iferr(hyperellcharpoly(Mod(1,p)*f), E, 0);
  if(chi == 0 || poldegree(chi) != 4 || !polisirreducible(chi), next);
  c3 = polcoef(chi,3); c2 = polcoef(chi,2);
  D = core(c3^2 - 4*(c2-2*p));
  if(D != 1, Ds = setunion(Ds, Set([D]))));
 print(cnames[i], " realfield-disc set: ", Ds));
}
quit
