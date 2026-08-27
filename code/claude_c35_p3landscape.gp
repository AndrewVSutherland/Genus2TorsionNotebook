\\ Step 1 for C35 target: (a) ALL genus-2 curves over F_3: can #J(F_3)=35 ?
\\ (b) Elkies A_1(5) chart residues at p=3,5,7: counts of 5|J, 7|J, 35|J.
default(parisize, 256000000);

\\ (a) exhaustive over F_3: every genus-2 curve /F_3 has a model y^2=f, deg f in {5,6}, f squarefree
{
p = 3;
cnt = vector(60);
tot = 0; sq = 0;
for(d = 5, 6,
  for(lc = 1, p-1,
    for(n = 0, p^d - 1,
      my(cs = digits(n, p));
      cs = concat(vector(d - #cs), cs);
      my(f = lc*'x^d + Pol(cs, 'x));
      tot++;
      my(fp = f * Mod(1,p));
      if (!issquarefree(fp), next);
      sq++;
      my(L = hyperellcharpoly(fp));
      my(J = subst(L, variable(L), 1));
      cnt[J]++;
    )
  )
);
print("F_3 exhaustive: models ", tot, "  squarefree ", sq);
print("attained #J(F_3) values and counts:");
for(J = 1, 60, if (cnt[J] > 0, print("  #J=", J, "  count ", cnt[J])));
print("#J=35 count: ", cnt[35]);
print("multiples of 5 attained: ", [J | J <- [5,10,15,20,25,30,35,40,45,50,55], cnt[J]>0]);
print("multiples of 7 attained: ", [J | J <- [7,14,21,28,35,42,49,56], cnt[J]>0]);
}

\\ (b) Elkies chart over F_p, p = 3, 5, 7
chart(p) = {
  my(sm=0, c5=0, c7=0, c35=0, bad=0, tot=0, x='x);
  for(a0=0,p-1, for(a1=0,p-1, for(a2=0,p-1,
    tot++;
    my(Q = Mod(1,p)*(a2*x^2+a1*x+a0), Qp = Q - x, H = Qp - x*Q, fz = H^2 + 4*Q^2*Qp);
    if (poldegree(fz) < 5 || !issquarefree(fz), bad++; next);
    my(L = hyperellcharpoly(fz), J = subst(L, variable(L), 1));
    sm++;
    if (J%5==0, c5++);
    if (J%7==0, c7++);
    if (J%35==0, c35++; print("   35-hit p=",p," triple ",[a0,a1,a2]," #J=",J));
  )));
  print("chart p=",p,": total ",tot," nongenus2 ",bad," smooth ",sm," | 5|J: ",c5," 7|J: ",c7," 35|J: ",c35);
}
chart(3); chart(5); chart(7);
quit;
