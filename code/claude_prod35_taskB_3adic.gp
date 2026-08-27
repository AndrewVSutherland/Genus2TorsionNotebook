\\ Task B: 3-adic sharpening for the Z/35 A_1(5) sieve.
\\ Part 1: exact genus2red verdict at p=3 for every GENUINE raw survivor of H=64
\\         (measures what a sharpened 3-adic filter would kill on its own).
\\ Part 2: valuation-pattern classifier: which v3(q0,q1,q2) patterns admit
\\         good reduction at 3 (necessarily with the magic Weil polynomial).
default(parisize, 512000000);
MAGIC = x^4 + 4*x^3 + 9*x^2 + 12*x + 9;   \\ unique F_3 Weil poly with 35 | #J

mkf(q0,q1,q2) = {
  my(Q = q2*'x^2 + q1*'x + q0, Qp = Q - 'x, H = Qp - 'x*Q);
  H^2 + 4*Q^2*Qp;
}

\\ verdict at 3: [class, detail]
\\ class 1 = C good reduction, magic poly (PASS - rigid)
\\ class 2 = C good reduction, wrong poly (KILL)
\\ class 3 = J bad reduction at 3 (conservative pass; detail = NU type string)
\\ class 4 = J good (e3=0) but C bad: compact type E1xE2 (needs traces {-1,-3})
verd3(fz) = {
  my(G = genus2red(fz, 3));
  my(FaN = G[2], PQ = G[3], V = G[4]);
  my(e3 = 0);
  for(i=1, matsize(FaN)[1], if(FaN[i,1]==3, e3 = FaN[i,2]));
  if (e3 == 0,
    my(P = PQ[1]*Mod(1,3), Q = PQ[2]*Mod(1,3));
    my(L = iferr(hyperellcharpoly([P,Q]), err, 0));
    if (L == 0, return([4, V[3][1]]));
    if (L == MAGIC*Mod(1,3)*1 || L == MAGIC, return([1, Str(L)]), return([2, Str(L)]));
  );
  return([3, V[3][1]]);
}

{
\\ ---------- Part 1 ----------
S = readstr("/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/prod35/h64_all.txt");
ngen = 0; nc = vector(4); types = Map();
for(i = 1, #S,
  my(w = strsplit(S[i], " "));
  if (#w < 4, next);
  my(q0 = eval(w[2]), q1 = eval(w[3]), q2 = eval(w[4]));
  my(fz = mkf(q0,q1,q2));
  if (poldegree(fz) < 5 || !issquarefree(fz), next);
  ngen++;
  my(v = verd3(fz));
  nc[v[1]]++;
  if (v[1] >= 2 && v[1] != 3,
    my(k = concat(Str(v[1]), concat(":", v[2])));
    mapput(types, k, if(mapisdefined(types,k), mapget(types,k), 0)+1));
);
print("PART1: genuine survivors H=64: ", ngen);
print("  class1 Cgood+magic: ", nc[1], "  class2 Cgood+wrong(KILL): ", nc[2], "  class3 Jbad: ", nc[3], "  class4 Jgood-Cbad(E1xE2): ", nc[4]);
my(M = Mat(types)); for(i=1, matsize(M)[1], print("   ", M[i,1], "  x", M[i,2]));
\\ ---------- Part 2 ----------
print("PART2: valuation pattern classifier, v in {-2..2}^3, 40 samples each");
NS = 40;
for(v0 = -2, 2, for(v1 = -2, 2, for(v2 = -2, 2,
  my(good = 0, magic = 0, bad = 0, degen = 0);
  for(s = 1, NS,
    my(u0 = 0, u1 = 0, u2 = 0);
    while(u0 % 3 == 0, u0 = random(200)-100); while(u1 % 3 == 0, u1 = random(200)-100); while(u2 % 3 == 0, u2 = random(200)-100);
    my(q0 = u0*3^v0, q1 = u1*3^v1, q2 = u2*3^v2);
    my(fz = mkf(q0,q1,q2));
    if (poldegree(fz) < 5 || !issquarefree(fz), degen++; next);
    my(vv = verd3(fz));
    if (vv[1] == 1, good++; magic++);
    if (vv[1] == 2, good++);
    if (vv[1] == 3, bad++);
  );
  if (good > 0,
    printf("  v=(%d,%d,%d): good=%d (magic=%d) bad=%d degen=%d\n", v0,v1,v2, good, magic, bad, degen));
)));
print("PART2 done (patterns with good==0 omitted)");
}
quit
