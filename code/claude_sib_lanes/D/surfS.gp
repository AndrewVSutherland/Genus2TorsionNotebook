\\ surfS.gp — Lane D: direct search of the T5 hit locus via the surface S.
\\ For member rho' = rn/rd, a (2,2,2,12) hit at u REQUIRES V'4 square, i.e.
\\ z^2 = F(u) := (q rho' - 1)(q rho' - 2u + 1), q = 4u^2-6u+3  — the curve C_rho'
\\ (genus 1).  On C_rho', V'4 == V'1 V'2 V'3 mod squares, so a point of C_rho'
\\ with X1, X2, X3 squares is a FULL HIT.  We enumerate C_rho' points by height
\\ (hyperellratpoints) and exact-test X1..X4 — reaches u-heights far beyond sweeps.
\\ Integer-cleared conditions (from t5sweep.c, u=p/q):
X123(p, q, rn, rd) = {
  my(QQ = 4*p^2-6*p*q+3*q^2,
     C3 = 16*p^4-40*p^3*q+40*p^2*q^2-18*p*q^3+3*q^4,
     C2 = -16*p^4+32*p^3*q-28*p^2*q^2+10*p*q^3-q^4,
     C1 = (8*p^3-12*p^2*q+10*p*q^2-3*q^3)*q,
     C0 = (-2*p+q)*q^3,
     X1, X2, X3, X4);
  X1 = rn*(QQ*rn-(2*p-q)*q*rd);
  X2 = QQ*rn^2-(4*p^2-4*p*q+2*q^2)*rn*rd+(2*p-q)*q*rd^2;
  X3 = rn*(C3*rn^3+C2*rn^2*rd+C1*rn*rd^2+C0*rd^3);
  X4 = C3*rn^2+(-16*p^3+28*p^2*q-18*p*q^2+4*q^3)*q*rn*rd+(2*p-q)^2*q^2*rd^2;
  [X1, X2, X3, X4];
}
sqpos(x) = (x > 0) && issquare(x);
\\ P(u) integer-cleared model of C_rho' : y^2 = (q*rn - rd)(q*rn - (2u-1)*rd)
Pmem(rn, rd) = my(u='x, q=4*u^2-6*u+3); (q*rn - rd)*(q*rn - (2*u-1)*rd);
\\ verbose scan (known members)
scanmember(rn, rd, H) = {
  my(P = Pmem(rn, rd), pts, np=0, nhit=0);
  if (poldisc(P) == 0, print("member ", rn, "/", rd, ": degenerate P, skip"); return);
  pts = hyperellratpoints(P, H);
  for (i = 1, #pts,
    my(u = pts[i][1], p = numerator(u), q = denominator(u), t);
    if (p == 0 || p == q || 2*p == q, next);
    t = X123(p, q, rn, rd);
    if (t[3] == 0 || t[4] == 0, next);   \\ degenerate section points
    np++;
    if (sqpos(t[1]) && sqpos(t[2]) && sqpos(t[3]),
      nhit++;
      print("HIT-CANDIDATE member ", rn, "/", rd, " u=", p, "/", q,
            "  X4sq=", sqpos(t[4])));
  );
  print("member ", rn, "/", rd, " H=", H, ": C-points(nondeg) = ", np, ", full = ", nhit);
}
\\ quiet scan for batch runs: prints only hits; returns [#nondeg pts, #hits]
qscan(rn, rd, H) = {
  my(P = Pmem(rn, rd), pts, np=0, nhit=0);
  if (poldisc(P) == 0, return([-1, 0]));
  pts = hyperellratpoints(P, H);
  for (i = 1, #pts,
    my(u = pts[i][1], p = numerator(u), q = denominator(u), t);
    if (p == 0 || p == q || 2*p == q, next);
    t = X123(p, q, rn, rd);
    if (t[3] == 0 || t[4] == 0, next);
    np++;
    if (sqpos(t[1]) && sqpos(t[2]) && sqpos(t[3]),
      nhit++;
      print("HIT-CANDIDATE member ", rn, "/", rd, " u=", p, "/", q,
            "  X4sq=", sqpos(t[4]))));
  [np, nhit];
}
