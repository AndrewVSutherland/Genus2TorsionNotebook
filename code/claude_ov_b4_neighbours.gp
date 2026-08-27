/* claude_ov_b4_neighbours.gp -- dump every Flynn(6,2) member and its
   Galois-stable Richelot neighbour(s) as integral models, for the Magma
   exact-torsion stage.  Output: results/claude_ov_b4_neighbours.txt with
   lines   n|sg|seedtype|seed2rank|nbtype|nb2rank|guard|seedpoly|nbpoly       */

read("/home/claude/torsion_jac/code/claude_ov_b4_richelot.gp");
if(type(NLIM) != "t_INT", NLIM = 40);
if(type(ROWS) != "t_STR", ROWS = "/home/claude/torsion_jac/results/claude_ov_b4_flynn62_rows.txt");
if(type(OUT)  != "t_STR", OUT  = "/home/claude/torsion_jac/results/claude_ov_b4_neighbours.txt");

Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;

{
  my(lines, nseed=0, nnb=0, nraise=0, stats=Map());
  system(Str("rm -f ", OUT));
  lines = Vecrev(readstr(ROWS));   /* smallest |n| first */
  foreach(lines, ln,
    my(rw, n, sg, t0, f, bl, fdegs, ftr);
    if(#ln < 3, next);
    rw = eval(ln);
    n = rw[1]; sg = rw[2]; t0 = rw[5];
    if(abs(n) > NLIM, next);
    f = Fl(t0);
    if(poldisc(f) == 0, next);
    nseed++;
    fdegs = b4_factortype(f); ftr = b4_tworank(fdegs);
    bl = b4_blocks(f);
    foreach(bl, br,
      my(rr = b4_resroots(br[2]));
      for(j = 1, #rr,
        my(o = b4_richelot(f, br[1], br[2], j), g, gd, gt);
        if(o[1] != 1,
           print("DEGEN n=", n, " sg=", sg, " j=", j, " code=", o[1]); next);
        g = b4_intmodel(o[2]);
        gd = b4_factortype(g); gt = b4_tworank(gd);
        nnb++;
        mapput(stats, Str(gd), if(mapisdefined(stats,Str(gd)), mapget(stats,Str(gd)), 0)+1);
        if(gt > ftr, nraise++;
           print("RAISE n=", n, " sg=", sg, " seed2rank=", ftr, " nb2rank=", gt,
                 " nbtype=", gd));
        write(OUT, Str(n, "|", sg, "|", fdegs, "|", ftr, "|", gd, "|", gt, "|",
                       o[5], "|", b4_intmodel(f), "|", g));
        print("NB n=", n, " sg=", sg, " seedtype=", fdegs, " nbtype=", gd,
              " nb2rank=", gt, " guard=", o[5],
              " nbdigits=", #Str(vecmax(apply(abs, Vec(g))))))));
  print("STATS neighbour factor types:");
  foreach(Vec(Set(Vec(Mat(stats)[,1]))), k, print("   ", k, " -> ", mapget(stats,k)));
  print("SEEDS ", nseed, "  NEIGHBOURS ", nnb, "  RAISES ", nraise);
  print("SEARCH_DONE b4neighbours NLIM=", NLIM);
}
quit
