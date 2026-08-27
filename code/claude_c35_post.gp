\\ Post-process height-8 survivors: drop globally degenerate triples,
\\ then extend necessary-condition sieve (35 | #J(F_q) at good q) through q <= 97.
default(parisize, 256000000);
S = readstr("c35_search_h16.txt");
PL = [31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
{
degen = 0; genuine = 0; killed = vector(#PL); finals = List();
for(i = 1, #S,
  my(w = strsplit(S[i], " "));
  if (#w < 4, next);
  my(q0 = eval(w[2]), q1 = eval(w[3]), q2 = eval(w[4]));
  my(Q = q2*'x^2 + q1*'x + q0, Qp = Q - 'x, H = Qp - 'x*Q);
  my(fz = H^2 + 4*Q^2*Qp);
  if (poldegree(fz) < 5 || !issquarefree(fz), degen++; next);
  genuine++;
  my(den = 1); for(k=0, poldegree(fz), den = lcm(den, denominator(polcoef(fz,k))));
  my(alive = 1);
  for(j = 1, #PL,
    my(q = PL[j]);
    if (den % q == 0, next);
    my(fq = fz * Mod(1,q));
    if (poldegree(fq) < 5 || !issquarefree(fq), next);  \\ bad reduction: conservative pass
    my(L = hyperellcharpoly(fq), J = subst(L, variable(L), 1));
    if (J % 35 != 0, killed[j]++; alive = 0; break);
  );
  if (alive, listput(finals, [q0,q1,q2]));
);
print("degenerate ", degen, "  genuine ", genuine);
print("kills by prime: ", [strjoin([Str(PL[j]),":",Str(killed[j])]) | j <- [1..#PL]]);
print("final survivors: ", #finals);
for(i = 1, #finals, print("FINAL ", finals[i]));
}
quit;
