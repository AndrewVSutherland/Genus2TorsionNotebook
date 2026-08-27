\\ Lane 8 (overnight 2026-07-25): isogeny-class ("Weil polynomial") landscape for the
\\ HLP split targets, over a long range of primes.
\\
\\ For each prime p and each target group G, enumerate every integral Weil polynomial
\\     chi(T) = T^4 - s1 T^3 + s2 T^2 - p s1 T + p^2      (all roots of modulus sqrt p)
\\ compatible with  G  contained in  J(Q)_tors  for a curve with good reduction at p:
\\     #G | #J(F_p) = p^2 + 1 - (p+1)s1 + s2      and    (T-1)^r | chi mod l  for (Z/l)^r in G.
\\ Report how many of those classes are ABSOLUTELY SIMPLE (root-power / D4 criterion).
\\
\\ Interpretation.  This is an upper bound on what can occur (not every Weil class is a
\\ Jacobian -- Howe-Nart-Ritzenthaler); the exhaustive CURVE enumeration in
\\ claude_ov_lane8_landscape.c settles realizability for p <= 23.  A prime with ZERO
\\ compatible classes forces bad reduction there; a prime with compatible classes but no
\\ absolutely simple one can never serve as a certifying prime.
\\
\\ usage:  gp -q code/claude_ov_lane8_weilscan.gp

abssimple(chi) = {
  my(v);
  if(!polisirreducible(chi), return(0));
  for(n = 2, 12,
    v = polresultant(subst(chi, x, y), x - y^n, y);
    v = v / polcoeff(v, poldegree(v));
    if(poldegree(v) != 4 || !polisirreducible(v), return(0));
  );
  1;
}

isweil(chi, p) = {
  my(r = polroots(chi), s = sqrt(p));
  for(i = 1, #r, if(abs(abs(r[i]) - s) > 1e-9*s, return(0)));
  1;
}

{TARGETS = [
  ["Z/35",     35, [[5,1],[7,1]]],
  ["Z/45",     45, [[3,1],[5,1]]],
  ["Z/7xZ/7",  49, [[7,2]]],
  ["Z/5xZ/10", 50, [[5,2],[2,1]]],
  ["Z/63",     63, [[7,1],[3,1]]],
  ["Z/70",     70, [[2,1],[5,1],[7,1]]],
  ["Z/2xZ/24", 48, [[2,2],[3,1]]],
  ["Z/3xZ/12", 36, [[3,2],[2,1]]],
  ["Z/8xZ/8",  64, [[2,2]]],
  ["Z/2xZ/4xZ/8", 64, [[2,3]]],
  ["Z/40 CTRL", 40, [[2,1],[5,1]]],
  ["Z/30 CTRL", 30, [[2,1],[3,1],[5,1]]],
  ["Z/2xZ/2xZ/12 CTRL", 48, [[2,3],[3,1]]]
];}

rankok(chi, l, r) = if(r == 0, 1, ((chi * Mod(1,l)) % ((x-1)^r)) == 0);

{
  my(PMAX = 200, res = matrix(#TARGETS, 60), plist = List(), npr = 0);
  forprime(p = 3, PMAX,
    npr++;
    listput(plist, p);
    my(S = ceil(4*sqrt(p)));
    for(t = 1, #TARGETS,
      my(N = TARGETS[t][2], reqs = TARGETS[t][3], ncomp = 0, nsimp = 0, exmp = 0);
      for(s1 = -S, S,
        for(s2 = -2*p - 1, s1^2\4 + 2*p + 1,
          my(nJ = p^2 + 1 - (p+1)*s1 + s2);
          if(nJ % N != 0 || nJ <= 0, next);
          my(chi = x^4 - s1*x^3 + s2*x^2 - p*s1*x + p^2);
          if(!isweil(chi, p), next);
          my(good = 1);
          for(j = 1, #reqs,
            if(p == reqs[j][1], next);
            if(!rankok(chi, reqs[j][1], reqs[j][2]), good = 0; break));
          if(!good, next);
          ncomp++;
          if(abssimple(chi), nsimp++; if(exmp == 0, exmp = [s1, s2, nJ]));
        )
      );
      res[t, npr] = [ncomp, nsimp, exmp];
    )
  );
  printf("Weil-polynomial (isogeny class) landscape, p = 3 .. %d\n", PMAX);
  printf("For each target: 'c/s' = (#compatible classes)/(#absolutely simple among them)\n\n");
  printf("%-20s", "prime");
  for(k = 1, npr, printf("%8d", plist[k])); printf("\n");
  for(t = 1, #TARGETS,
    printf("%-20s", TARGETS[t][1]);
    for(k = 1, npr, printf("%8s", Str(res[t,k][1], "/", res[t,k][2])));
    printf("\n");
  );
  printf("\nFIRST prime admitting an absolutely simple compatible class, and that class:\n");
  for(t = 1, #TARGETS,
    my(done = 0);
    for(k = 1, npr,
      if(res[t,k][2] > 0 && done == 0,
        done = 1;
        printf("  %-20s p = %3d   (s1,s2) = (%d,%d)  #J = %d\n",
               TARGETS[t][1], plist[k], res[t,k][3][1], res[t,k][3][2], res[t,k][3][3]));
    );
    if(done == 0, printf("  %-20s NONE below %d  <<<<<< LOCAL SIMPLICITY OBSTRUCTION\n", TARGETS[t][1], PMAX));
  );
  printf("\nPrimes with ZERO compatible classes (forced bad reduction):\n");
  for(t = 1, #TARGETS,
    my(bad = List());
    for(k = 1, npr, if(res[t,k][1] == 0, listput(bad, plist[k])));
    printf("  %-20s %s\n", TARGETS[t][1], Str(Vec(bad)));
  );
  quit;
}
