\\ Post-process sieve HIT lines (file name in the global variable HITFILE):
\\ exact CF order of D_inf over Q, factor type, discriminant.
\\ Lines: "HIT w a b c d".
x='x;
sqrtpolypart(f)={my(s=x^3,dd); for(k=1,3, dd=f-s^2; if(poldegree(dd)<=2,break); s=s+(polcoeff(dd,6-k)/(2*polcoeff(s,3)))*x^(3-k)); s;}
cforder(f,ms)={my(s,Pi,Qi,total,ai,Pn,Qn); s=sqrtpolypart(f); Pi=0*x; Qi=1+0*x; total=0;
  for(i=0,ms, if(Qi==0,return(0)); ai=(Pi+s)\Qi; total+=poldegree(ai); Pn=ai*Qi-Pi;
    if((f-Pn^2)%Qi!=0, return(0)); Qn=(f-Pn^2)\Qi; Pi=Pn; Qi=Qn;
    if(i>=1 && poldegree(Qi)<=0 && Qi!=0, return(total))); 0;}
chartf(v) = x*(x-v[1])*(x^2+v[2]*x+v[3])*(x^2+v[4]*x+v[5]);
{
  data = readstr(HITFILE);
  nh = 0; nord11 = 0;
  for(k=1,#data,
    pieces = strsplit(data[k], " ");
    if(#pieces != 6, next);
    if(pieces[1] != "HIT", next);
    v = vector(5, j, eval(pieces[j+1]));
    f = chartf(v);
    dd = poldisc(f);
    if(dd == 0, next);
    nh++;
    ord = cforder(f,40);
    fa = factor(f); tp = vecsort(vector(matsize(fa)[1],j,poldegree(fa[j,1])));
    if(ord == 11, nord11++);
    print("CAND ", v, " CFord=", ord, " type=", tp, " disc=", factor(dd)));
  print("CHECK_DONE lines=", nh, " cford11=", nord11);
}
