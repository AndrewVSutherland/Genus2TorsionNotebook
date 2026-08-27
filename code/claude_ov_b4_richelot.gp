/* claude_ov_b4_richelot.gp -- Lane 4 (route B4): Galois-stable Richelot walker.

   Given a sextic f = c * F1 * F2 with F1 a MONIC rational quadratic block and
   F2 the complementary MONIC quartic whose Galois group lies in D4 (resolvent
   cubic has a rational root), the partition {c*F1, G2, G3} of the six
   Weierstrass points -- where F2 = G2*G3 over a quadratic field -- is
   GALOIS-STABLE, so the (2,2)-kernel and the Richelot isogeny are defined
   over Q even though the kernel is not pointwise rational.

   Richelot: f = G1 G2 G3, Delta = det(coeff matrix),
             H_i = G_j' G_k - G_j G_k',  codomain  Delta y^2 = H1 H2 H3,
             i.e.  y^2 = (H1 H2 H3)/Delta.   (Delta and H1H2H3 are each
             ANTI-invariant under the conjugation swapping G2,G3, so the
             quotient is rational; the product Delta*H1H2H3 is the quadratic
             twist by A and gives WRONG torsion.)

   All arithmetic in Q[w]/(w^2 - A) (t_POLMOD); the output sextic is rational.
   Quartic splitting used:  F2 depressed = X^4 + p X^2 + q X + r, resolvent
   z^3 - p z^2 - 4 r z + (4pr - q^2); for a rational root z0 put A = z0 - p:
     A != 0 :  (X^2 + aX + b)(X^2 - aX + c), a^2 = A,
               c = ((p+A) + q/a)/2, b = ((p+A) - q/a)/2;
     A == 0 :  forces q = 0; (X^2 + b)(X^2 + c), b+c = p, bc = r, i.e.
               b,c = (p -+ sqrt(p^2-4r))/2  (the biquadratic case).

   Self-test at the bottom reproduces the Magma-verified 2-rank-1 control.
*/

x='x; w='w; z='z;

b4_depress(q4) = {
  my(a3,a2c,a1,a0,s,p,q,r);
  a3 = polcoef(q4,3,x); a2c = polcoef(q4,2,x);
  a1 = polcoef(q4,1,x); a0 = polcoef(q4,0,x);
  s = a3/4;
  p = a2c - 3*a3^2/8;
  q = a1 - a3*a2c/2 + a3^3/8;
  r = a0 - a1*a3/4 + a2c*a3^2/16 - 3*a3^4/256;
  [s,p,q,r];
};

b4_resolvent(q4) = {
  my(d, p, q, r);
  d = b4_depress(q4); p = d[2]; q = d[3]; r = d[4];
  z^3 - p*z^2 - 4*r*z + (4*p*r - q^2);
};

b4_resroots(q4) = {
  my(c, f, out=List());
  c = b4_resolvent(q4);
  if(poldegree(c) < 1, return([]));
  f = factor(c);
  for(i=1, matsize(f)[1],
    if(poldegree(f[i,1]) == 1,
      listput(out, -polcoef(f[i,1],0,z)/polcoef(f[i,1],1,z))));
  Vec(out);
};

b4_factortype(f) = {
  my(fa);
  fa = factor(f);
  vecsort(apply(pp->poldegree(pp), fa[,1]~));
};

/* 2-rank of J(Q) for y^2 = f, deg f = 6 squarefree, from the factor degrees */
b4_tworank(degs) = {
  my(k, neven);
  /* deg-5 models carry the extra rational Weierstrass point at infinity */
  if(vecsum(Vec(degs)) == 5, degs = vecsort(concat(Vec(degs), [1])));
  k = #degs;
  neven = 0;
  for(m = 0, 2^k - 1,
    my(s = 0);
    for(i = 1, k, if(bittest(m, i-1), s += degs[i]));
    if(s % 2 == 0, neven++));
  round(log(neven)/log(2)) - 1;
};

/* All rational quadratic BLOCKS of f: [monic quadratic B, monic quartic R]
   with f = pollead(f) * B * R.  Sources: an irreducible rational quadratic
   factor, or a product of two rational linear factors. */
b4_blocks(f) = {
  my(fa, c, lins=List(), quads=List(), out=List(), R, mf, chk);
  if(poldegree(f) != 6, return([]));
  fa = factor(f);
  c = pollead(f);
  /* RECORDED TRAP: PARI factor() returns PRIMITIVE INTEGRAL factors and drops
     the content, so the factors are NOT monic and their product is NOT f.
     Re-monicise every factor and verify f = pollead(f) * prod(monic factors). */
  chk = c;
  for(i=1, matsize(fa)[1],
    if(fa[i,2] > 1, return([]));          /* not squarefree */
    mf = fa[i,1]/pollead(fa[i,1]);
    chk *= mf^fa[i,2];
    if(poldegree(mf) == 1, listput(lins, mf));
    if(poldegree(mf) == 2, listput(quads, mf)));
  if(chk != f, error("b4_blocks: monic reconstruction failed"));
  foreach(Vec(quads), B,
    R = f/(c*B); if(type(R)=="t_POL" && poldegree(R)==4 && pollead(R)==1, listput(out, [B,R])));
  for(i=1, #lins-1, for(j=i+1, #lins,
    my(B = lins[i]*lins[j]);
    R = f/(c*B); if(type(R)=="t_POL" && poldegree(R)==4 && pollead(R)==1, listput(out, [B,R]))));
  Vec(out);
};

/* Richelot step.  f: rational sextic; B: monic rational quadratic block;
   R: complementary monic quartic; j: index of the rational resolvent root.
   Returns [1, fprime, A, Delta] or [0] (degenerate) or [-1/-2, data] (bug). */
b4_richelot(f, B, R, j) = {
  my(rts, d, s, p, q, r, A, A2, MW, ww, bb, cc, G1, G2, G3, chk,
     m, Delta, H1, H2, H3, fp, fpl, c);
  c = pollead(f);
  if(f - c*B*R != 0, return([-3]));
  rts = b4_resroots(R);
  if(#rts < j, return([0]));
  d = b4_depress(R); s = d[1]; p = d[2]; q = d[3]; r = d[4];
  A = rts[j] - p;
  if(A != 0,
    MW = w^2 - A;
    ww = Mod(w, MW);
    cc = ((p + A) + q*ww/A)/2;
    bb = ((p + A) - q*ww/A)/2;
    G2 = (x+s)^2 + ww*(x+s) + bb;
    G3 = (x+s)^2 - ww*(x+s) + cc;
  ,
    if(q != 0, return([0]));
    A2 = p^2 - 4*r;
    if(A2 == 0, return([0]));
    MW = w^2 - A2;
    ww = Mod(w, MW);
    bb = (p - ww)/2; cc = (p + ww)/2;
    G2 = (x+s)^2 + bb;
    G3 = (x+s)^2 + cc;
  );
  G1 = Mod(1,MW)*(c*B);
  chk = G1*G2*G3 - Mod(1,MW)*f;
  if(chk != 0, return([-1, lift(chk)]));
  m = [polcoef(G1,2,x), polcoef(G1,1,x), polcoef(G1,0,x);
       polcoef(G2,2,x), polcoef(G2,1,x), polcoef(G2,0,x);
       polcoef(G3,2,x), polcoef(G3,1,x), polcoef(G3,0,x)];
  Delta = matdet(m);
  if(Delta == 0, return([0]));
  H1 = deriv(G2,x)*G3 - G2*deriv(G3,x);
  H2 = deriv(G3,x)*G1 - G3*deriv(G1,x);
  H3 = deriv(G1,x)*G2 - G1*deriv(G2,x);
  /* Codomain: Delta*y^2 = H1H2H3 up to the cyclic-order convention for H_i.
     Two traps here, BOTH measured against Magma's RichelotIsogenousSurfaces
     and against isogeny-invariance of #J(F_p):
       (a) Delta*H1H2H3 is the quadratic TWIST BY A of the right answer (A is
           a non-square exactly when the kernel is Galois-stable but not
           rational), so one must DIVIDE by Delta, not multiply;
       (b) the sign convention H_i = G_j'G_k - G_jG_k' vs G_k'G_j - G_kG_j'
           flips the product (three factors), i.e. twists by -1.
     The guard below re-derives the twist from #J(F_p) rather than trusting
     the convention. */
  fp = -H1*H2*H3/Delta;
  fpl = lift(fp);
  if(poldegree(fpl, w) > 0, return([-2, fpl]));
  fpl = subst(fpl, w, 0);
  if(poldegree(fpl) < 5 || poldisc(fpl) == 0, return([0]));
  /* twist guard: a Richelot isogeny preserves #J(F_p) at every good prime. */
  my(good = 0, okp = 0, okm = 0);
  foreach([11,13,17,19,23,29,31,37,41,43,47,53], pp,
    if(good >= 4, next);
    my(Ns = b4_jcount(f,pp), Np = b4_jcount(fpl,pp), Nm = b4_jcount(-fpl,pp));
    if(Ns == 0 || Np == 0 || Nm == 0, next);
    good++;
    if(Ns == Np, okp++);
    if(Ns == Nm, okm++));
  if(good == 0, return([1, fpl, lift(A), lift(Delta), 0]));
  if(okp == good && okm < good, return([1, fpl, lift(A), lift(Delta), good]));
  if(okm == good && okp < good, return([1, -fpl, lift(A), lift(Delta), good]));
  if(okp == good && okm == good, return([1, fpl, lift(A), lift(Delta), -good]));
  return([-4, fpl, okp, okm, good]);
};

/* y^2 = f  ->  an EQUIVALENT integral model.  f = c*g with g integral
   primitive and c = a/b in lowest terms; then (b*y)^2 = a*b*g, so the model is
   a*b*g.  NO integer factorization is used (factoring the content of a
   thousand-digit model is what made the first version hang); square factors
   are only stripped by trial division over small primes to keep sizes down. */
b4_intmodel(f) = {
  my(c, g, a, b, m, e, q);
  if(f == 0, return(f));
  c = content(f);
  g = f/c;
  a = numerator(c); b = denominator(c);
  m = a*b;
  forprime(pp = 2, 5000,
    e = valuation(m, pp);
    if(e >= 2, m /= pp^(2*(e\2))));
  m * g;
};

/* #J(F_p) for y^2 = f (odd good p), via hyperellcharpoly */
b4_jcount(f, p) = {
  my(g, D);
  if(p == 2, return(0));
  g = b4_intmodel(f);
  D = poldisc(g);
  if(D == 0 || D % p == 0, return(0));
  if(Mod(pollead(g),p) == 0 && poldegree(g) == 6, return(0));
  subst(hyperellcharpoly(Mod(1,p)*g), x, 1);
};

/* squarefree cores of c3^2-4(c2-2p) over the first np good primes (RM screen) */
b4_discsig(f, np) = {
  my(S=List(), p=3, n=0, cp, c3, c2, d, D, g);
  g = b4_intmodel(f);
  D = poldisc(g);
  while(n < np && p < 600,
    if(D % p != 0 && Mod(pollead(g),p) != 0,
      cp = hyperellcharpoly(Mod(1,p)*g);
      c3 = polcoef(cp,3); c2 = polcoef(cp,2);
      d = c3^2 - 4*(c2 - 2*p);
      if(d != 0, listput(S, core(d)); n++));
    p = nextprime(p+1));
  Set(Vec(S));
};

b4_selftest() = {
  my(f, bl, out, degs, ok);
  /* Magma-verified control: TwoTorsionSubgroup = [2], ONE Richelot neighbour
     over Q from a non-split Galois-stable kernel (claude_richelot_rank1.m) */
  f = (x^2+x+1)*(x^4-4*x^2+2);
  print("SELFTEST seed type=", b4_factortype(f), " 2-rank=", b4_tworank(b4_factortype(f)));
  bl = b4_blocks(f);
  print("SELFTEST blocks: ", #bl);
  ok = 0;
  foreach(bl, br,
    my(rr = b4_resroots(br[2]));
    print("SELFTEST  block=", br[1], "  resolvent rational roots=", rr);
    for(j=1, #rr,
      out = b4_richelot(f, br[1], br[2], j);
      if(out[1] == 1,
        degs = b4_factortype(out[2]);
        print("SELFTEST  -> codomain f'=", b4_intmodel(out[2]));
        print("SELFTEST     type=", degs, " 2-rank=", b4_tworank(degs),
              "  A=", out[3], " Delta=", out[4]);
        ok = 1
      , print("SELFTEST  -> degenerate ", out))));
  /* a second control with a genuinely rational (2,2)-split (2-rank 2) */
  f = (x^2+x+1)*(x^2+2*x+2)*(x^2+3*x+5);
  bl = b4_blocks(f);
  foreach(bl, br,
    my(rr = b4_resroots(br[2]));
    for(j=1, #rr,
      out = b4_richelot(f, br[1], br[2], j);
      print("SELFTEST2 rational-split seed 2-rank=", b4_tworank(b4_factortype(f)),
            " block=", br[1], " -> ",
            if(out[1]==1, Str(b4_intmodel(out[2])), Str(out)))));
  ok;
};
