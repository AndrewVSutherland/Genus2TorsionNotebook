SetColumns(0);
if not assigned q then q := 31; elif Type(q) eq MonStgElt then q := StringToInteger(q); end if;
kq := GF(q);
Fx<X> := RationalFunctionField(kq);
Py<Y> := PolynomialRing(Fx);
Q8 := Y^8 + (216*X^4+72*X^3-24*X^2)*Y^4
    + (-1296*X^6-1728*X^5-432*X^4+64*X^3)*Y^2
    + (-3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4);
F<yy> := FunctionField(Q8);
xx := F!X;
yox := yy/xx;
pl1 := Places(F,1);
printf "q=%o  genus %o  #pl1 %o\n", q, Genus(F), #pl1;
for i in [1..#pl1] do
  R := pl1[i];
  vx := Valuation(xx, R);
  if vx lt 0 then
    cy := Evaluate(yox, R);
    printf "%3o INF  vx=%o  y/x=%o  v(y/x-c)=%o  v(1/x)=%o\n", i, vx, cy, Valuation(yox-F!cy,R), Valuation(1/xx,R);
  else
    cx := Evaluate(xx, R); cy := Evaluate(yy, R);
    printf "%3o FIN  x=%o v(x-c)=%o  y=%o v(y-c)=%o\n", i, cx, Valuation(xx-F!cx,R), cy, Valuation(yy-F!cy,R);
  end if;
end for;
print "SIGDIAG_DONE";
quit;
