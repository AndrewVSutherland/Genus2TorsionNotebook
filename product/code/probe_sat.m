SetColumns(0);
SetClassGroupBounds("GRH");
E := EllipticCurve([0,0,1,-1,0]);   // 37a, rank 1, generator (0,0)
G := E![0,0];
S := [2*G];
t0 := Realtime();
T := Saturation(S, 20);
printf "Saturation([2G],20): %o pts, %o s, ht(first)=%o (ht(G)=%o, ht(2G)=%o)\n",
    #T, Realtime()-t0, Height(T[1]), Height(G), Height(2*G);
S2 := [2*G];
T2 := Saturation(S2, 3);
printf "Saturation([2G],3): ht(first)=%o\n", Height(T2[1]);
printf "PROBE_SAT_DONE\n";
quit;
