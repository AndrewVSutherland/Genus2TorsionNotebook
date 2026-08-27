/* claude_ov_88_ds1_sieve.gp -- Lane 5 (2026-07-25).
   Enumerate the DOUBLE-STAGE-1 locus of Nicholls' Lambda_334 family in the new
   (m,n,p) chart derived this session:

     t     = (m^2+1)/(2m),  alpha = (m^2-1)/(2m)          [t^2-1 = alpha^2]
     s     = (n^2+alpha^4)/(2n)                            [s^2-alpha^4 = square]
     beta  = t(1-p^2)/(1+p^2), gamma = 2tp/(1+p^2)         [t^2-beta^2 = gamma^2]
     family/eta condition:  w^2 = s^2(1+p^2)^2 - 4 alpha^2 t^2 p^2
                              = (s p^2 + 2 alpha t p + s)(s p^2 - 2 alpha t p + s)

   The first two lines are exactly "class {0,a} in 2J2(Q)", the third is exactly
   "class {c,oo} in 2J2(Q)" (both verified against exact Magma torsion, 235/235),
   and the fourth is the condition that the point lies on the family at all.
   So every hit of this sieve is a DOUBLE-STAGE-1 member: both prescribed
   2-torsion classes of J2 are 2-divisible over Q -- the necessary condition for
   [8,8] on J1.  Previously such members could only be produced by a
   MordellWeilGroup computation on E_{s,t}: w^2 = z(z+s^2)(z+A).

   Integral form.  m = m1/m2, n = n1/n2, p = p1/p2, all reduced, all positive
   (only t^2, s^2, p^2 enter the curve, so signs are redundant), and
     M = 2 m1 m2, T = m1^2+m2^2, Al = m1^2-m2^2,
     Sn = n1^2 M^4 + Al^4 n2^2,   k0 = 4 Al T n1 n2 M^2,
     V  = Sn^2 (p1^2+p2^2)^2 - k0^2 p1^2 p2^2     must be a perfect square.

   Usage: gp -q -D timer=0 ... with MMAX,NMAX,PMAX,PART,NPARTS in the environment.
*/
{
MMAX = eval(getenv("MMAX"));
NMAX = eval(getenv("NMAX"));
PMAX = eval(getenv("PMAX"));
PART = eval(getenv("PART"));
NPARTS = eval(getenv("NPARTS"));
if(MMAX==0, MMAX=12); if(NMAX==0, NMAX=12); if(PMAX==0, PMAX=300);
if(NPARTS==0, NPARTS=1);

mlist = List();
for(m1=1, MMAX, for(m2=1, MMAX,
  if(gcd(m1,m2)==1 && m1!=m2, listput(mlist, [m1,m2]))));
nlist = List();
for(n1=1, NMAX, for(n2=1, NMAX,
  if(gcd(n1,n2)==1, listput(nlist, [n1,n2]))));
plist = List();
for(p2=2, PMAX, for(p1=1, p2-1,
  if(gcd(p1,p2)==1, listput(plist, [p1,p2]))));

nhit = 0; ntest = 0;
for(im=1, #mlist,
  if((im-1) % NPARTS != PART, next);
  m1 = mlist[im][1]; m2 = mlist[im][2];
  M = 2*m1*m2; T = m1^2+m2^2; Al = m1^2-m2^2;
  M4 = M^4; Al4 = Al^4; M2 = M^2;
  for(inn=1, #nlist,
    n1 = nlist[inn][1]; n2 = nlist[inn][2];
    Sn = n1^2*M4 + Al4*n2^2;
    k0 = 4*Al*T*n1*n2*M2;
    Sn2 = Sn^2; k02 = k0^2;
    for(ip=1, #plist,
      p1 = plist[ip][1]; p2 = plist[ip][2];
      pp = p1^2+p2^2;
      q = pp^2;
      \\ drop the spurious component gamma^2 = alpha^2 (<=> beta^2 = 1 <=> b = 1,
      \\ a DEGENERATE member: the quintic acquires a double root at x=1).
      if(4*T^2*p1^2*p2^2 == Al^2*q, next);
      V = Sn2*q - k02*p1^2*p2^2;
      ntest++;
      if(V > 0 && issquare(V),
        nhit++;
        print("DS1HIT ", m1, " ", m2, " ", n1, " ", n2, " ", p1, " ", p2)
      )
    )
  )
);
print("SIEVE_DONE part=", PART, " tested=", ntest, " hits=", nhit);
quit;
}
