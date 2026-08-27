/* claude_ov_erratum_control.gp -- genus2red ERRATUM propagation (2026-07-25).
 *
 * The PARI half of the control.  Same 14 LMFDB curves as
 * code/claude_ov_erratum_control.m, with their PUBLISHED conductors (the first
 * component of the LMFDB label).  We show that genus2red returns exactly the
 * ODD PART of the published conductor in every case, and that it always
 * reports the sentinel exponent -1 at p = 2.
 *
 * Usage:   gp -q code/claude_ov_erratum_control.gp
 * Markers: GPCTRL / GPCTRLSUMMARY / ERRATUM_GPCTRL_DONE
 */
{
ctrl = [
  ["295.a.295.2",     x^5-40*x^3+22*x^2+389*x-608,          x^2+x+1,  295],
  ["464.a.464.1",     -x^6-2*x^5-2*x^4-x^3,                 x+1,      464],
  ["464.a.29696.2",   4*x^5+33*x^4+72*x^3+16*x^2+x,         x,        464],
  ["1416.b.135936.1", -2*x^4-x^3+x+1,                       x^3+x,    1416],
  ["2600.a.338000.1", 10*x^5+8*x^4-5*x^3-3*x^2+x,           x,        2600],
  ["1575.c.1",        -27*x^5+23*x^4+8*x^3-4*x^2-x,         x^2+x,    1575],
  ["3942.b.3",        9*x^5+20*x^4-3*x^3-6*x^2-x,           x^2+x,    3942],
  ["9450.b.2",        -225*x^5+146*x^4-16*x^3-6*x^2+x,      x^2+x,    9450],
  ["10512.n.1",       -12*x^5+6*x^4+12*x^3-7*x^2,           x^2+1,    10512],
  ["12300.e.2",       x^6-7*x^5+8*x^4+16*x^3-5*x^2-5*x,     x^2+x,    12300],
  ["19044.h.2",       x^6-3*x^5+9*x^4-5*x^3+12*x^2-6*x,     x^2+x,    19044],
  ["39600.dq.1",      -6*x^5-14*x^4+18*x^3+11*x^2-12*x+2,   x^2+1,    39600],
  ["8136.c.1",        9*x^6-3*x^5-12*x^4+3*x^3+2*x^2,       x^2+1,    8136],
  ["28200.e.1",       9*x^6-3*x^5-42*x^4+34*x^3+4*x^2-3*x,  x^2+1,    28200]
];
nodd = 0; nsent = 0; ntot = #ctrl;
for(i=1, ntot,
  lab = ctrl[i][1]; f = ctrl[i][2]; h = ctrl[i][3]; Npub = ctrl[i][4];
  G = genus2red([f,h]);
  Ngp = G[1]; M = G[2];
  oddpub = Npub / 2^valuation(Npub,2);
  ok = (Ngp == oddpub);
  /* sentinel: does the factorisation matrix carry a row [2, -1]? */
  sent = 0;
  for(j=1, matsize(M)[1], if(M[j,1]==2 && M[j,2]==-1, sent = 1));
  if(ok, nodd++);
  if(sent, nsent++);
  print("GPCTRL ", lab, "  Npub=", Npub, "  v2(Npub)=", valuation(Npub,2),
        "  oddpart(Npub)=", oddpub, "  genus2red_N=", Ngp,
        "  EQUALS_ODDPART=", ok, "  sentinel[2,-1]=", sent);
);
print("GPCTRLSUMMARY  genus2red_N == oddpart(Npub) in ", nodd, " of ", ntot,
      " ; sentinel [2,-1] present in ", nsent, " of ", ntot);
print("ERRATUM_GPCTRL_DONE");
}
quit;
