/* Lane 8 (overnight 2026-07-25): postprocessor for claude_ov_lane8_landscape.c
 *
 * Reads a "s1 s2 count" histogram of the Frobenius data of EVERY genus-2 curve
 * y^2 = f(x) over F_p (deg f in {5,6}, leading coefficient in {1,nu}) and reports,
 * for each target torsion order N, the isogeny classes with N | #J(F_p), together
 * with whether each is simple / absolutely simple.
 *
 * Reading:  a curve C/Q whose Jacobian has a rational point of order N and which
 * has GOOD REDUCTION at p must reduce into one of the listed classes.  If none of
 * them is absolutely simple, p can never be a certifying prime for such a C.
 *
 * usage:
 *   gp -q code/claude_ov_lane8_landscape.gp \
 *      <(echo 'analyse(5, readhist("results/claude_ov_lane8_landscape_p5.hist"), [21,27,35,40,45,49,50,63,70]);')
 */

\\ absolutely-simple test (root-power / "D4" criterion): chi irreducible and the
\\ charpoly of alpha^n irreducible of degree 4 for n = 2..12.
abssimple(chi) = {
  my(v);
  if(!polisirreducible(chi), return(0));
  for(n=2, 12,
    v = polresultant(subst(chi,x,y), x - y^n, y);
    v = v/polcoeff(v, poldegree(v));
    if(poldegree(v) != 4 || !polisirreducible(v), return(0));
  );
  1;
}

readhist(fn) = {
  my(h=List(), s, v, l);
  s = readstr(fn);
  for(i=1, #s,
    l = s[i];
    if(l == "", next);
    if(Vecsmall(l)[1] == 35, next);
    v = eval(Str("[", strjoin(strsplit(l, " "), ","), "]"));
    if(#v == 3, listput(h, v));
  );
  Vec(h);
}

analyse(p, hist, targets) = {
  my(tot=0, nJ, chi, s1, s2, cnt);
  for(i=1, #hist, tot += hist[i][3]);
  printf("\n================ p = %d ================\n", p);
  printf("curves enumerated (deg f in {5,6}, lc in {1,nu}) = %d ; distinct isogeny classes = %d\n", tot, #hist);
  printf("Weil range for #J(F_p): [%.2f, %.2f]\n", (sqrt(p)-1)^4, (sqrt(p)+1)^4);
  for(t=1, #targets,
    my(N=targets[t], hitc=0, hitn=0, simpc=0, simpn=0, irrc=0, irrn=0, lines=List());
    for(i=1, #hist,
      s1 = hist[i][1]; s2 = hist[i][2]; cnt = hist[i][3];
      nJ = p^2 + 1 - (p+1)*s1 + s2;
      if(nJ % N == 0,
        chi = x^4 - s1*x^3 + s2*x^2 - p*s1*x + p^2;
        my(irr = polisirreducible(chi), as = abssimple(chi));
        hitc += cnt; hitn += 1;
        if(irr, irrc += cnt; irrn += 1);
        if(as, simpc += cnt; simpn += 1);
        listput(lines, [s1, s2, nJ, cnt, irr, as]);
      );
    );
    printf("\n--- N = %d : %d isogeny classes / %d curves (%.5f%%) ; chi irreducible: %d cls / %d curves ; ABSOLUTELY SIMPLE: %d cls / %d curves\n",
           N, hitn, hitc, 100.0*hitc/tot, irrn, irrc, simpn, simpc);
    for(i=1, #lines,
      my(L=lines[i]);
      printf("    s1=%4d s2=%5d  #J=%7d  curves=%9d  irred=%d  abs_simple=%d\n",
             L[1], L[2], L[3], L[4], L[5], L[6]);
    );
  );
}
