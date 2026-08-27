\\ probeB.gp - [2,22] leads on the two order-11 families.
\\ (1) generic factorization of Flynn F_t and Daowsud-Schmidt G_t over Q(t);
\\ (2) rational-quadratic-factor incidence: eliminate t from the two remainder
\\     coefficients of F_t mod (x^2+u*x+v); factor the (u,v)-projection W(u,v);
\\     for factors quadratic in u or v print the discriminant (hyperelliptic
\\     double-cover model of that component);
\\ (3) small (u,v) sweep via gcd of the two t-quadratics;
\\ (4) witness diagnostics: factor types + CF order of D_infinity for the three
\\     production generic witnesses 1192.a [22], 1312.c [22], 1416.b [2,14],
\\     and CF self-tests (f14->14, f18->18, Flynn->11, DS->11).
default(parisizemax, 2*10^9);
x = 'x; t = 't; u = 'u; v = 'v;

Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;
DS(tt) = x^6-4*x^5+8*(1+tt)*x^4-(10+32*tt)*x^3+8*(1+6*tt+2*tt^2)*x^2-4*(1+6*tt+16*tt^2)*x+64*tt^2+1;

print("Flynn factorization over Q(t):");
print(lift(factor(Fl(t))));
print("DS factorization over Q(t):");
print(lift(factor(DS(t))));

\\ ---------- CF order machinery over Q ----------
SqrtPolyPart(f) = {
 my(s = x^3);
 for(k=1, 3,
  my(d = f - s^2);
  if(poldegree(d) <= 2, break);
  s += polcoef(d, 6-k)/(2*polcoef(s,3)) * x^(3-k));
 s }
CFOrd(f, maxsteps) = {
 my(s = SqrtPolyPart(f), Pi = 0, Qi = 1, total = 0);
 for(i=0, maxsteps,
  if(Qi == 0, return(0));
  my(ai = divrem(Pi+s, Qi)[1]);
  total += poldegree(ai);
  my(Pn = ai*Qi - Pi);
  if((f - Pn^2) % Qi != 0, return(0));
  my(Qn = divrem(f - Pn^2, Qi)[1]);
  Pi = Pn; Qi = Qn;
  if(i >= 1 && poldegree(Qi) <= 0 && Qi != 0, return(total)));
 0 }

\\ self-tests
tv1 = CFOrd((x^2+1)*(x^4+5*x^2+4*x+4), 60);
tv2 = CFOrd((x^2-x+1)*(x^4-x^3+9*x^2+8*x-8), 60);
print("CF self-test: f14 -> ", tv1, " (want 14), f18 -> ", tv2, " (want 18)");
for(k=1, 3, my(tt=[2,5,-3][k]); print("Flynn t=",tt," CFOrd=", CFOrd(Fl(tt),60)));
for(k=1, 3, my(tt=[1,3,-2][k]); print("DS t=",tt," CFOrd=", CFOrd(DS(tt),60)));

\\ ---------- incidence analysis ----------
famnames = ["Flynn", "DS"];
{fams = [Fl(t), DS(t)];}
{for(fi=1, 2,
  my(f = fams[fi]);
  my(dv = divrem(f, x^2+u*x+v));
  my(r = dv[2]);
  my(e1 = polcoef(r,1,x), e0 = polcoef(r,0,x));
  my(W = polresultant(e1, e0, t));
  my(fw = factor(W));
  print("==== ", famnames[fi], ": t-resultant of quadratic-factor incidence ====");
  for(k=1, matsize(fw)[1],
    my(g = fw[k,1]);
    if(poldegree(g,u) == 0 && poldegree(g,v) == 0, next);
    print("  factor mult=", fw[k,2], " deg_u=", poldegree(g,u),
          " deg_v=", poldegree(g,v), " :");
    print("    ", g);
    if(poldegree(g,v) == 2,
      my(dsc = polcoef(g,1,v)^2 - 4*polcoef(g,2,v)*polcoef(g,0,v));
      print("    disc_v factored: ", factor(dsc)));
    if(poldegree(g,u) == 2,
      my(dsc = polcoef(g,1,u)^2 - 4*polcoef(g,2,u)*polcoef(g,0,u));
      print("    disc_u factored: ", factor(dsc)));
  );
  \\ (u,v) sweep
  my(HH = 10, found = 0);
  for(du=1, HH, for(nu=-HH*du, HH*du,
   if(gcd(nu,du) != 1, next);
   my(u0 = nu/du);
   for(dvv=1, HH, for(nv=-HH*dvv, HH*dvv,
    if(nv == 0 || gcd(nv,dvv) != 1, next);
    my(v0 = nv/dvv);
    my(p1 = subst(subst(e1,u,u0),v,v0), p0 = subst(subst(e0,u,u0),v,v0));
    if(p1 == 0 && p0 == 0, next);
    my(g = if(p1 == 0, p0, if(p0 == 0, p1, gcd(p1,p0))));
    if(type(g) != "t_POL" || poldegree(g,t) < 1, next);
    foreach(nfroots(,g), t0,
      my(ff = subst(f, t, t0));
      if(poldegree(ff,x) < 5 || poldisc(ff) == 0, next);
      found++;
      my(fa = factor(ff));
      my(degs = vecsort(apply(poldegree, fa[,1]~)));
      my(cford = if(issquare(pollead(ff)) && poldegree(ff)==6,
                    CFOrd(ff/pollead(ff), 60), -1));
      print("QF_POINT ", famnames[fi], " t=",t0," u=",u0," v=",v0,
            " type=",degs, " CFOrd=", cford);
    );
   ));
  ));
  print(famnames[fi], " sweep found=", found);
);}

\\ ---------- witness diagnostics ----------
{wit = [
 ["1192.a [22]",   4*(x^3-2*x^2-x+1) + (x^3+x)^2],
 ["1312.c [22]",   4*(x^6+4*x^5+7*x^4+5*x^3+2*x^2) + (x+1)^2],
 ["1416.b [2,14]", 4*(-2*x^4-x^3+x+1) + (x^3+x)^2],
 ["RM 19044.h [2,22]",     4*(x^6-3*x^5+9*x^4-5*x^3+12*x^2-6*x) + (x^2+x)^2],
 ["RM 152100.eb [2,2,14]", 4*(9*x^6+69*x^5+123*x^4-95*x^3-183*x^2+165*x-35) + (x^2+x)^2]
];}
{for(k=1, #wit,
  my(nm = wit[k][1], F = wit[k][2]);
  my(fa = factor(F));
  my(degs = vecsort(apply(poldegree, fa[,1]~)));
  my(lc = pollead(F));
  my(cford = if(poldegree(F)==6 && issquare(lc), CFOrd(F/lc, 60), -1));
  print("WITNESS ", nm, " type=", degs, " lc=", lc, " CFOrd(D_inf)=", cford);
);}
print("PROBE_B_DONE");
quit
