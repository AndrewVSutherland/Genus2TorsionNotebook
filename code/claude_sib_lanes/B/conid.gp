Ws(r) = {
  my(q = 4*'u^2-6*'u+3, u='u);
  [ r*(q*r-(2*u-1)),
    q*r^2-(4*u^2-4*u+2)*r+(2*u-1),
    r*((16*u^4-40*u^3+40*u^2-18*u+3)*r^3+(-16*u^4+32*u^3-28*u^2+10*u-1)*r^2+(8*u^3-12*u^2+10*u-3)*r+(-2*u+1)),
    (16*u^4-40*u^3+40*u^2-18*u+3)*r^2+(-16*u^3+28*u^2-18*u+4)*r+(4*u^2-4*u+1) ]
}
sqfreepart(p) = {my(f=factor(p), s=1); for(i=1,matsize(f)[1], if(f[i,2]%2==1, s*=f[i,1])); s}
{
forstep(i=1,6,1,
  my(r=[-49/240,289/240,-1/143,-25/551,-169/1431,841/697][i]);
  my(W=Ws(r));
  print("r = ",r);
  print("  sqfree(W1*W2*W4) = ", sqfreepart(W[1]*W[2]*W[4]));
  print("  sqfree(W1*W4)    = ", sqfreepart(W[1]*W[4]));
  print("  sqfree(W2*W4)    = ", sqfreepart(W[2]*W[4]));
  print("  sqfree(W1*W2)    = ", sqfreepart(W[1]*W[2]));
  print("  sqfree(W1*W2*W3) = ", sqfreepart(W[1]*W[2]*W[3]));
);
}
quit
