/* claude_ov_b4_c43tables.gp -- Lane 4 (route B4): mod-p Legendre tables for the
   ONLY 2-rank-raise locus of the (4,3) component of Flynn's quadratic-factor
   incidence variety.

   Setting (all derived and verified in code/claude_ov_b4_c43param.m):
     (4,3) component, clean parametrization by w:
        u(w) = (2w^3-2w^2+2w-1)/w^2,  v(w) = ((w^2-w+1)/w)^2,
        t(w) = -((w-1)(w^2-w+1)/w^2)^2,
     F_{t(w)} = L(x) * G1(x) * C(x)  with degrees 1,2,3 over Q(w)  => 2-rank 1.
     Writing T = (w^3-2w^2+2w-1)/w^2,  the cubic is
        C(x) = x^3 + (1-2T)x^2 + (1+T^2)x - T^2 = (x-1)T^2 - 2x^2 T + x(x^2+x+1),
     so C has a rational root iff x is a square, x = rho^2, and then
        T = rho(rho^2-rho+1)/(rho-1)   [the sign flip is rho -> -rho].
     With z = 1/(1-rho) this reads  T(w) + T(z) = 0, i.e. the raise locus is
        C3 : (w^3-2w^2+2w-1) z^2 + (z^3-2z^2+2z-1) w^2 = 0        (genus 3)
     which is SYMMETRIC in (w,z).  In p = w+z, q = wz it becomes
        p^2 - p(q^2+2q) + (4q^2-2q) = 0,
     a genus-1 curve  y^2 = q^4+4q^3-12q^2+8q  (y = 2p - q^2 - 2q), i.e.
        E : V^2 = X^3 - 12X^2 + 32X + 64,   X = 8/q,  V = 8y/q^2,
     which is 92.a1 AGAIN -- the same rank-1 conductor-92 curve that carries the
     (6,2) component.  E(Q) = Z.G with G = (8,8), trivial torsion.

   So: C3(Q)  <->  { n in Z : delta(nG) = p^2 - 4q  is a rational square }.
   This script tabulates, for each prime, which residues n mod ord(G mod p) are
   allowed by the Legendre symbol of delta.

   Output: results/claude_ov_b4_c43tab_<TAG>.txt, blocks "p N_p" + one mask line
   (format consumed by code/claude_ov_b4_sieve.c with nloci = 1).

   Usage: echo 'PMAX=3000;TAG="main";read("code/claude_ov_b4_c43tables.gp")' | gp -q
*/

RESDIR = if(type(RESDIR)=="t_STR", RESDIR, "/home/claude/torsion_jac/results");
PMAX = if(type(PMAX)=="t_INT", PMAX, 3000);
TAG  = if(type(TAG)=="t_STR", TAG, "main");

EA = [0,-12,0,32,64];
GX = 8; GY = 8;

/* exact rational delta at n*G (for the self-test) */
{
  my(E = ellinit(EA), G = [GX,GY]);
  print("E cond = ", ellglobalred(E)[1], "  torsion = ", elltors(E)[1]);
  for(n = -6, 6,
    my(P = ellmul(E,G,n), qq, yy, pp, dd);
    if(#P < 2, print("  n=", n, " : O  (q=0 branch, delta=0, DEGENERATE)"); next);
    if(P[1] == 0, print("  n=", n, " : X=0 -> q=oo  DEGENERATE"); next);
    qq = 8/P[1]; yy = 8*P[2]/P[1]^2; pp = ((qq^2+2*qq)+yy)/2; dd = pp^2-4*qq;
    print("  n=", n, " q=", qq, " y=", yy, " p=", pp, " delta=", dd,
          "  square? ", issquare(dd));
  );
}

/* ------------------------------------------------------------------ tables */
{
  my(fn, plist, nb = 0);
  fn = strprintf("%s/claude_ov_b4_c43tab_%s.txt", RESDIR, TAG);
  write(fn, "# claude_ov_b4_c43tables.gp  locus: delta = p^2-4q is a square on E=92.a1, G=(8,8)");
  plist = primes([5, PMAX]);
  foreach(plist, p,
    if(p == 23, next);                       /* bad reduction (cond 92 = 4*23) */
    my(e, G, N, P, s, qq, yy, pp, dd, ok);
    e = ellinit(EA, p);
    if(type(e) != "t_VEC" || #e == 0, next);
    G = [Mod(GX,p), Mod(GY,p)];
    if(!ellisoncurve(e, G), next);
    N = ellorder(e, G);
    if(N == 0 || N > 4*p, next);
    s = vector(N);
    P = [0];                                  /* identity */
    for(r = 0, N-1,
      ok = 1;
      if(#P >= 2,
        if(P[1] == 0,
          ok = 1                              /* q = oo : degenerate, allow */
        ,
          qq = 8/P[1]; yy = 8*P[2]/P[1]^2;
          pp = ((qq^2+2*qq)+yy)/2; dd = pp^2-4*qq;
          ok = if(kronecker(lift(dd), p) == -1, 0, 1);
        );
      );
      s[r+1] = ok;
      P = elladd(e, P, G);
    );
    write(fn, p, " ", N);
    write(fn, concat(apply(c -> if(c, "1", "0"), Vec(s))));
    nb++;
  );
  print("TABLES WRITTEN: ", nb, " primes, PMAX=", PMAX, " -> ", fn);
}
quit;
