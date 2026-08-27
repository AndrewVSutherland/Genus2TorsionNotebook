/* claude_ov_b4_walk.gp -- Lane 4 (route B4): Galois-stable Richelot BFS from
   the Flynn (6,2) [2,4]-type stream (infinitely many generic [22] seeds).

   Reads results/claude_ov_b4_flynn62_rows.txt (produced by
   claude_ov_b4_flynn62.gp), rebuilds each member's sextic F_t, and walks the
   Galois-stable Richelot graph to depth DEPTH, reporting every node whose
   2-rank is >= 2 (a [2,22] candidate) together with the odd-torsion gate
   (11 | #J(F_p)) and the RM disc-census signature.

   Usage:
     gp -q -D parisize=... -D nbthreads=1 <<EOF
       ROWS="/path/rows.txt"; DEPTH=3; NLIM=40; read(".../claude_ov_b4_walk.gp");
     EOF
*/

read("/home/claude/torsion_jac/code/claude_ov_b4_richelot.gp");

if(type(DEPTH) != "t_INT", DEPTH = 3);
if(type(NLIM)  != "t_INT", NLIM  = 40);
if(type(ROWS)  != "t_STR", ROWS = "/home/claude/torsion_jac/results/claude_ov_b4_flynn62_rows.txt");

Fl(tt) = x^6+2*x^5+(2*tt+3)*x^4+2*x^3+(tt^2+1)*x^2+2*tt*(1-tt)*x+tt^2;

/* one Richelot expansion: all Galois-stable (2,2)-kernels reachable via a
   rational quadratic block; returns a vector of codomain sextics */
b4_children(f) = {
  my(out=List(), bl, seenc=Map());
  bl = b4_blocks(f);
  foreach(bl, br,
    my(rr = b4_resroots(br[2]));
    for(j=1, #rr,
      my(o = b4_richelot(f, br[1], br[2], j));
      if(o[1] == 1,
        my(g = b4_intmodel(o[2]), k);
        k = Str(g);
        if(!mapisdefined(seenc, k), mapput(seenc, k, 1); listput(out, g)))));
  Vec(out);
};

{
  my(fh, lines, nseed, nnode, nhit, tot_children, stats, GLOBALHITS);
  lines = readstr(ROWS);
  print("ROWS ", #lines, "  DEPTH=", DEPTH, "  NLIM=", NLIM);
  nseed = 0; nnode = 0; nhit = 0; tot_children = 0;
  stats = Map();
  GLOBALHITS = List();
  foreach(lines, ln,
    my(rw, n, sg, uu, vv, t0, f, seen, frontier, d);
    if(#ln < 3, next);
    rw = eval(ln);
    n = rw[1]; sg = rw[2]; uu = rw[3]; vv = rw[4]; t0 = rw[5];
    if(abs(n) > NLIM, next);
    f = Fl(t0);
    if(poldisc(f) == 0, next);
    if(f % (x^2+uu*x+vv) != 0, print("BADSEED n=",n); next);
    nseed++;
    /* seed sanity: type [2,4] or [1,1,4], 2-rank 1, 11 | #J(F_p) */
    my(dg = b4_factortype(f), tr = b4_tworank(dg), okodd = 1);
    foreach([13,17,19,29,31], pp,
      my(N = b4_jcount(f,pp)); if(N != 0 && N % 11 != 0, okodd = 0));
    if(tr != 1 || okodd != 1,
      print("SEEDWARN n=",n," sg=",sg," type=",dg," 2rank=",tr," odd11=",okodd));
    seen = Map(); mapput(seen, Str(b4_intmodel(f)), 1);
    frontier = [f];
    for(d = 1, DEPTH,
      my(nxt = List());
      foreach(frontier, g,
        my(ch = b4_children(g));
        tot_children += #ch;
        foreach(ch, h,
          my(k = Str(h), degs, trk, N11);
          if(mapisdefined(seen, k), next);
          mapput(seen, k, 1);
          nnode++;
          degs = b4_factortype(h); trk = b4_tworank(degs);
          mapput(stats, Str(degs), if(mapisdefined(stats,Str(degs)), mapget(stats,Str(degs)), 0) + 1);
          if(trk >= 2,
            nhit++;
            N11 = [b4_jcount(h,pp) : pp <- [13,17,19,23,29,31,37]];
            print("HIT seed_n=", n, " sg=", sg, " depth=", d, " type=", degs,
                  " 2rank=", trk, " Jcounts=", N11);
            print("HITPOLY ", h);
            listput(GLOBALHITS, [n, sg, d, h]);
          );
          listput(nxt, h)));
      frontier = Vec(nxt);
      if(#frontier == 0, break));
  );
  print("STATS codomain factor types:");
  foreach(Vec(Set(Vec(Mat(stats)[,1]))), k, print("   ", k, " -> ", mapget(stats,k)));
  print("SEEDS ", nseed, "  NODES ", nnode, "  CHILD_EDGES ", tot_children,
        "  HITS ", nhit);
  if(#GLOBALHITS > 0,
    system("rm -f /home/claude/torsion_jac/results/claude_ov_b4_walk_hits.txt");
    foreach(Vec(GLOBALHITS), h,
      write("/home/claude/torsion_jac/results/claude_ov_b4_walk_hits.txt", h)));
  print("SEARCH_DONE b4walk DEPTH=", DEPTH, " NLIM=", NLIM);
}
quit
