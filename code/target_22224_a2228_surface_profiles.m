//////////////////////////////////////////////////////////////////////
// Point-count residue profiles on the genuine rational surface
//
//   (a,b,c,d)=(1,r,s,rs),
//   q=2t/(1-t^2),  s=(1-q^2*r)/(q^2-r).
//
// For fixed t=m/n, A=(n^2-m^2)^2 and B=(2mn)^2, a projective integral
// tuple is (B-Ar, r(B-Ar), A-Br, r(A-Br)).  Each row records whether a
// specialization r in P^1(F_p) has boundary reduction or 3 | #J(F_p).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned PrimeList then
    PrimeList := [37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,
                  107,109,113,127,131,137,139,149,151,157,163,167,173,
                  179,181,191,193,197,199];
elif Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s):s in Split(PrimeList,",")];
end if;
if not assigned TList then TList := [<1,2>,<2,3>]; end if;
if not assigned output_file then
    output_file := "results/target_22224_a2228_surface_profiles.tsv";
end if;
if not assigned log_file then
    log_file := "results/target_22224_a2228_surface_profiles.log";
end if;

SetLogFile(log_file : Overwrite := true);
Z := Integers(); out := Open(output_file,"w");
fprintf out,"t_num\tt_den\tp\tr\tboundary\tjorder\tjmod3\tallowed\n";

for tr in TList do
    tm,tn := Explode(tr);
    for p in PrimeList do
        F := GF(p); P<x> := PolynomialRing(F);
        m:=F!tm; n:=F!tn;
        A := (n^2-m^2)^2; B := (2*m*n)^2;
        allowed_count:=0; boundary_count:=0; good_allowed:=0;
        for ri in [0..p] do
            if ri eq p then
                // r=infinity: leading coefficients in r.
                vals := [-A,-A,-B,-B]; label := "inf";
            else
                r:=F!ri; sd:=B-A*r; sn:=A-B*r;
                vals := [sd,r*sd,sn,r*sn]; label:=IntegerToString(ri);
            end if;
            allzero := &and[z eq 0:z in vals];
            boundary := allzero or &or[z eq 0:z in vals]
                        or #Set([z^2:z in vals]) ne 4;
            nord:=0; jmod:=-1;
            if boundary then
                boundary_count +:= 1;
            else
                f:=x*&*[x+z^2:z in vals];
                assert Degree(f) eq 5 and Discriminant(f) ne 0;
                nord:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
                jmod:=nord mod 3;
            end if;
            allowed:=boundary or jmod eq 0;
            if allowed then allowed_count+:=1; end if;
            if not boundary and allowed then good_allowed+:=1; end if;
            fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                tm,tn,p,label,boundary select 1 else 0,nord,jmod,
                allowed select 1 else 0;
        end for;
        print "SURFACE_PROFILE","t",tr,"p",p,"boundary",boundary_count,
              "allowed",allowed_count,"good_allowed",good_allowed;
    end for;
end for;

delete out;
print "TARGET_22224_A2228_SURFACE_PROFILES_DONE",output_file;
UnsetLogFile();
quit;
