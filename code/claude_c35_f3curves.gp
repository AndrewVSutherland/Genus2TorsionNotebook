\\ List the 12 genus-2 models over F_3 with #J(F_3)=35 and factor their L-polys.
{
p = 3;
for(d = 5, 6,
  for(lc = 1, p-1,
    for(n = 0, p^d - 1,
      my(cs = digits(n, p));
      cs = concat(vector(d - #cs), cs);
      my(f = lc*'x^d + Pol(cs, 'x));
      my(fp = f * Mod(1,p));
      if (!issquarefree(fp), next);
      my(L = hyperellcharpoly(fp));
      my(J = subst(L, variable(L), 1));
      if (J == 35,
        my(fa = factor(L));
        print("f = ", lift(fp), "  charpoly = ", lift(L), "  factors ", matsize(fa)[1],
              "  irr4: ", if(matsize(fa)[1]==1 && poldegree(fa[1,1])==4, "YES", "no"));
      );
    )
  )
);
}
quit;
