\\ claude_ov_88obs.gp -- Lane 5 (2026-07-25, third session).
\\
\\ THE LOCAL READING OF THE [8,8] FAILURE ON Lambda_334.
\\
\\ On the double-stage-1 locus the sextic has factor type [1,1,2,2], so
\\     L = Q[x]/f1 = Q x Q x K1 x K3,   K1 = Q(sqrt(D1)),  K3 = Q(sqrt(a)).
\\ An order-4 class D is 2-divisible in J(Q) iff there is a single c in Q^* with
\\ u_D(theta_i) in c*(K_i^*)^2 for every component.  Writing u_D = lambda_i * square
\\ in K_i, and using  Q^* cap (K_i^*)^2 = (Q^*)^2 u d_i (Q^*)^2,  the two RATIONAL
\\ components force c = lambda_lin exactly, so the test is
\\     def1 := lambda_lin*lambda_1 in {1, D1}   and   def3 := lambda_lin*lambda_3 in {1, a}
\\ modulo squares.  This script measures the STRICTLY WEAKER, purely local relaxation
\\     def1 in N(K1^*)  and  def3 in N(K3^*)
\\ i.e. solubility of the conics  z^2 - D1*w^2 = def1*y^2  and  z^2 - a*w^2 = def3*y^2,
\\ which by Hasse-Minkowski is everywhere-local solubility, and reports the exact set
\\ of obstructing places when it fails.  If even this relaxation fails, the failure of
\\ the [8,8] route on the DS1 locus is LOCAL at the level of the lambda-compatibility
\\ conic, not a subtle global (Sha-like) failure.
\\
\\ Input : a file of "D1 a def1 def3" records (squarefree integers), one per line,
\\         default the records extracted from results/claude_ov_88verify.log.
\\ Usage : gp -q -D recfile=<path> code/claude_ov_88obs.gp

if(type(recfile) != "t_STR", recfile = "/home/claude/torsion_jac/data/claude_ov_88obs_recs.txt");

\\ obstructing places for "x is a norm from Q(sqrt(d))"; -1 denotes the real place
obs(x,d) = {
  my(P=List(), S=Set(concat([2], concat(factor(abs(x))[,1]~, factor(abs(d))[,1]~))));
  if(x<0 && d<0, listput(P,-1));
  for(i=1,#S, my(p=S[i]); if(hilbert(x,d,p)==-1, listput(P,p)));
  Vec(P);
}

{
  my(f=fileopen(recfile), ln, tot=0, ex1=0, ex3=0, nrm1=0, nrm3=0, nrmboth=0);
  my(c1=Map(), c3=Map(), hist1=Map(), minp=List(), soluble=List());
  while((ln=filereadstr(f))!=0,
    if(#ln==0 || Vecsmall(ln)[1]==35, next);
    my(w = apply(eval, Vec(strsplit(ln," "))));
    my(D1=w[1], a=w[2], d1=w[3], d3=w[4]);
    tot++;
    if(d1==1 || d1==D1, ex1++);
    if(d3==1 || d3==a, ex3++);
    my(o1=obs(d1,D1), o3=obs(d3,a));
    if(#o1==0, nrm1++);
    if(#o3==0, nrm3++);
    if(#o1==0 && #o3==0, nrmboth++);
    if(#o1==0 || #o3==0, listput(soluble, [D1,a,d1,d3,#o1,#o3]));
    for(i=1,#o1, my(p=o1[i]); mapput(c1,p, if(mapisdefined(c1,p),mapget(c1,p),0)+1));
    for(i=1,#o3, my(p=o3[i]); mapput(c3,p, if(mapisdefined(c3,p),mapget(c3,p),0)+1));
    mapput(hist1,#o1, if(mapisdefined(hist1,#o1),mapget(hist1,#o1),0)+1);
    if(#o1>0, listput(minp, vecmin([if(o1[i]==-1,0,o1[i]) | i<-[1..#o1]])));
  );
  fileclose(f);
  print("OBS_RECORDS ", tot);
  print("EXACT_GATE_PASS comp1=", ex1, " comp3=", ex3);
  print("LOCALLY_SOLUBLE comp1=", nrm1, " comp3=", nrm3, " both=", nrmboth);
  print("HIST_obstructing_places_comp1 ", [[k, mapget(hist1,k)] | k <- vecsort(Vec(hist1))]);
  my(ks=vecsort(Vec(c1)));
  print("NDISTINCT_obstruction_primes_comp1 ", #ks);
  print("TOP_obstruction_primes_comp1 ", vecsort([[mapget(c1,ks[i]), ks[i]] | i<-[1..#ks]],,4)[1..min(15,#ks)]);
  my(ks3=vecsort(Vec(c3)));
  print("NDISTINCT_obstruction_primes_comp3 ", #ks3);
  print("TOP_obstruction_primes_comp3 ", vecsort([[mapget(c3,ks3[i]), ks3[i]] | i<-[1..#ks3]],,4)[1..min(15,#ks3)]);
  if(#minp, print("SMALLEST_obstructing_place_comp1 max over records = ", vecmax(Vec(minp)),
                  "  (0 = the real place)"));
  print("PARTIALLY_SOLUBLE_RECORDS [D1,a,def1,def3,#obs1,#obs3]:");
  for(i=1,#soluble, print("  ", soluble[i]));
  print("OBS_DONE");
}
quit
