\\ Independent PARI validation for claude_z31_box31.c (Task B5).
\\ Run from the repo root:  gp -q code/claude_z31_box31_check.gp
\\ Reads results/claude_z31_box31_chkdata.gp (written by `claude_z31_box31 selftest`),
\\ recomputes N1 = #C(F_p), N2 = #C(F_p^2), J = #J(F_p) via hyperellcharpoly,
\\ and brute-force point counts (for the first curves) to pin the conventions.

read("results/claude_z31_box31_chkdata.gp");

\\ brute-force #C(F_p): affine points + 2*issquare(c6) at infinity
bruteN1(f, p) = {
  my(c = if(issquare(Mod(polcoeff(f, 6), p)), 2, 0));
  for(x = 0, p - 1,
    my(v = Mod(subst(f, variable(f), x), p));
    c += if(v == 0, 1, if(issquare(v), 2, 0)));
  c;
}
\\ brute-force #C(F_p^2) (curves in chkdata have p !| c6, so infinity gives 2)
bruteN2(f, p) = {
  my(w = ffgen([p, 2], 'w), c = 2);
  for(a = 0, p - 1, for(b = 0, p - 1,
    my(v = subst(f, variable(f), a + b*w));
    c += if(v == 0, 1, if(issquare(v), 2, 0))));
  c;
}

npass = 0; nfail = 0;
{for(i = 1, #chkdata,
  my(v = chkdata[i], f, p, P, J, a, b, N1, N2);
  f = Pol(v[1..7]); p = v[8];
  P = hyperellcharpoly(Mod(1, p)*f);
  J = subst(P, variable(P), 1);
  a = -polcoeff(P, 3); b = polcoeff(P, 2);
  N1 = p + 1 - a;
  N2 = p^2 + 1 - (a^2 - 2*b);
  if(N1 == v[9] && N2 == v[10] && J == v[11], npass++,
     nfail++; print("GPFAIL ", v, " gp:", [N1, N2, J]));
)}
print("GPCHECK charpoly pass=", npass, " fail=", nfail);

\\ brute-force convention pin: first 3 curves, N1 at all 6 primes, N2 at p=37,41
bpass = 0; bfail = 0;
{for(i = 1, min(18, #chkdata),
  my(v = chkdata[i], f, p, n1);
  f = Pol(v[1..7]); p = v[8];
  n1 = bruteN1(f, p);
  if(n1 == v[9], bpass++, bfail++; print("BRUTEFAIL N1 ", v, " brute:", n1));
  if(p <= 41,
    my(n2 = bruteN2(f, p));
    if(n2 == v[10], bpass++, bfail++; print("BRUTEFAIL N2 ", v, " brute:", n2)));
)}
print("GPCHECK brute pass=", bpass, " fail=", bfail);

\\ the RM witness: expect 31 | #J(F_p) at every good prime in the list
F = -3356*x^6 + 11364*x^5 - 18347*x^4 + 17202*x^3 - 9863*x^2 + 3264*x - 504;
print("WITNESS disc = ", factor(poldisc(F)));
{foreach([37, 41, 43, 47, 53, 59], p,
  if(Mod(poldisc(F), p) == 0 || Mod(polcoeff(F, 6), p) == 0,
     print("WITNESS p=", p, " BAD"),
     my(P = hyperellcharpoly(Mod(1, p)*F), J = subst(P, variable(P), 1));
     print("WITNESS p=", p, " J=", J, " Jmod31=", J % 31)));
}
quit;
