x='x;
{
N = -3744*x^4 + 117536*x^3 - 1300544*x^2 + 5379456*x - 3432704;
Dn = 625*x^4 - 15000*x^3 + 135000*x^2 - 540000*x + 810000;
S = N*Dn;
ker = 117*x^4 - 3673*x^3 + 40642*x^2 - 168108*x + 107272;
q = subst(S,x,0)/subst(ker,x,0);
c0 = core(numerator(q)*denominator(q));
print("twist constant core: ", c0);
model = c0*ker;
print("model: Y^2 = ", model);
pts = hyperellratpoints(model, 200000);
print("points H<=2*10^5: ", pts);
\\ if genus-1 with a point: rank
if(#pts > 0,
  P0 = pts[1];
  E = ellfromeqn('y^2 - subst(model,x,'x));
  print("ellfromeqn: ", E[1..5]);
  Ee = ellinit(E);
  print("cond = ", ellglobalred(Ee)[1], "  analytic rank = ", ellanalyticrank(Ee)[1], "  torsion = ", elltors(Ee)[1]));
}
quit
