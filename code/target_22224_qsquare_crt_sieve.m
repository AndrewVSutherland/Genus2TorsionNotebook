//////////////////////////////////////////////////////////////////////
// Second-stage finite enumeration for [2,2,2,24].
//
// Exhaustively enumerates the affine q=(x+t)^2 slice over each prime,
// with both rho and sigma sheets and all three D4 pairings.  Every open
// halving point is retained, checked in the finite Jacobian, and assigned
// the relative (M,N)-Jacobian rank used for Hensel lifting.
//
// Run from the repository root:
//   magma -b code/target_22224_qsquare_crt_sieve.m
//
// Options:
//   PrimeList:=29,31,37,41
//   output_stem:=results/target_22224_qsquare_crt_sieve
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned PrimeList then
    PrimeList := [29,31,37,41];
elif Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s):s in Split(PrimeList,",")];
end if;

if not assigned output_stem then
    output_stem := "results/target_22224_qsquare_crt_sieve";
end if;
if not assigned log_file then
    log_file := output_stem cat ".log";
end if;

SetLogFile(log_file : Overwrite:=true);

Z:=Integers();

function MakeMonic(g)
    return g/LeadingCoefficient(g);
end function;

function DirectData(K,aa,bb,cc,dd,uu,tt,vv)
    P<X>:=PolynomialRing(K);
    f:=X*(X+aa^2)*(X+bb^2)*(X+cc^2)*(X+dd^2);
    q:=X^2+uu*X+tt^2;
    h:=X^3+(1+3*uu)/2*X^2+vv*X+tt^3;
    return P,f,q,h;
end function;

// Return the marked D12 together with the CRT auxiliaries needed by the
// incidence cover.  The boolean means the stated affine CRT chart is open.
function MarkedD12Aux(K,aa,bb,cc,dd,uu,tt,vv)
    P,f,q,h:=DirectData(K,aa,bb,cc,dd,uu,tt,vv);
    X:=P.1;
    if h^2-f-q^3 ne 0 then
        return false,P!0,P!0,P!0,P!0,P!0,P!0,P!0,
               K!0,K!0,K!0,P!0;
    end if;

    g4:=(X-aa*bb)*(X-cc*dd)-X*(aa+bb)*(cc+dd);
    L4:=(X-aa*bb)*(cc+dd)+(aa+bb)*(X-cc*dd);
    if Degree(g4) ne 2 then
        return false,f,q,h,g4,L4,P!0,P!0,K!0,K!0,K!0,P!0;
    end if;
    v3:=h mod q;
    v4:=(-X*L4) mod g4;

    Rm:=(q*X) mod g4;
    Rn:=q mod g4;
    Rt:=(v4-v3) mod g4;
    a0:=Coefficient(Rm,0); a1:=Coefficient(Rm,1);
    b0:=Coefficient(Rn,0); b1:=Coefficient(Rn,1);
    c0:=Coefficient(Rt,0); c1:=Coefficient(Rt,1);
    det:=a0*b1-a1*b0;
    if det eq 0 then
        return false,f,q,h,g4,L4,P!0,P!0,K!0,K!0,det,P!0;
    end if;
    mc:=(c0*b1-c1*b0)/det;
    nc:=(a0*c1-a1*c0)/det;
    ell:=v3+q*(mc*X+nc);
    if mc eq 0 or (ell^2-f) mod (q*g4) ne 0 then
        return false,f,q,h,g4,L4,P!0,P!0,mc,nc,det,ell;
    end if;
    Uraw:=ExactQuotient(ell^2-f,q*g4);
    if Degree(Uraw) ne 2 then
        return false,f,q,h,g4,L4,P!0,P!0,mc,nc,det,ell;
    end if;
    U:=MakeMonic(Uraw);
    V:=(-ell) mod U;
    if (V^2-f) mod U ne 0 then
        return false,f,q,h,g4,L4,U,V,mc,nc,det,ell;
    end if;
    return true,f,q,h,g4,L4,U,V,mc,nc,det,ell;
end function;

function HalvingData(U,V,f)
    K:=BaseRing(Parent(f));
    A<M,N>:=PolynomialRing(K,2);
    AX<X>:=PolynomialRing(A);
    phi:=hom<Parent(f)->AX|X>;
    UX:=phi(U); VX:=phi(V); fX:=phi(f);
    ellH:=VX+UX*(M*X+N);
    if (ellH^2-fX) mod UX ne 0 then
        return false,A!0,A!0,AX!0,AX!0;
    end if;
    S:=ExactQuotient(ellH^2-fX,UX);
    s4:=Coefficient(S,4); s3:=Coefficient(S,3);
    s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
    E1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
    E0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
    return true,E1,E0,S,ellH;
end function;

function EvalCoeffPolynomial(g,vals,P)
    if IsZero(g) then return P!0; end if;
    return P![Evaluate(Coefficient(g,i),vals):i in [0..Degree(g)]];
end function;

function FInt(z)
    return Z!z;
end function;

function PolynomialKey(g,n)
    return Join([Sprint(FInt(Coefficient(g,i))):i in [0..n]],",");
end function;

function MarkedKey(f,U,V)
    return PolynomialKey(f,5) cat "|" cat PolynomialKey(U,2) cat
           "|" cat PolynomialKey(V,1);
end function;

function HalfKey(f,U,V,G,W)
    return MarkedKey(f,U,V) cat "|" cat PolynomialKey(G,2) cat
           "|" cat PolynomialKey(W,1);
end function;

function IsDivisibleBy2Finite(J,D)
    G,phi:=AbelianGroup(J);
    a:=D @@ phi;
    coords:=Eltseq(a);
    invs:=Invariants(G);
    for i in [1..#coords] do
        if GCD(2,invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

procedure Bump(~AA,key)
    if IsDefined(AA,key) then AA[key]+:=1; else AA[key]:=1; end if;
end procedure;

procedure EnumeratePrime(p,stem)
    K:=GF(p);
    points_name:=stem cat "_p" cat IntegerToString(p) cat "_points.tsv";
    boundary_name:=stem cat "_p" cat IntegerToString(p) cat "_boundary.tsv";
    mask_name:=stem cat "_p" cat IntegerToString(p) cat "_masks.tsv";
    pts:=Open(points_name,"w");
    bdy:=Open(boundary_name,"w");
    msk:=Open(mask_name,"w");
    fprintf pts,"p\tA\tB\tC\trho\tsigma\tpairing\ta\tb\tc\td\tmc\tnc\tU1\tU0\tV1\tV0\tM\tN\tG1\tG0\tW1\tW0\tordD\tordH\trelrankMN\n";
    fprintf bdy,"p\tA\tB\tC\trho\tsigma\tpairing\tM\tN\treason\n";
    fprintf msk,"p\tA\tB\tR\tS\trho_zero\tsigma_zero\tbranch_collision\tsmooth_slice_sheets\topen_presentations\ttarget_presentations\ttarget_halves\thalf_boundary_points\n";

    torus_pairs:=(p-1)^2;
    double_square_pairs:=0;
    rho_zero_pairs:=0; sigma_zero_pairs:=0; sheet_boundary_pairs:=0;
    affine_sheet_points:=0; collision_sheet_points:=0; smooth_sheet_points:=0;
    open_presentations:=0; chart_failures:=0; contact_failures:=0;
    order_tests:=0; order12_count:=0; order_other_count:=0;
    divisibility_tests:=0; finite_group_decompositions:=0;
    divisible_presentations:=0;
    target_presentations:=0; target_halves:=0; verified_halves:=0;
    half_boundary_points:=0; projective_MN_halves:=0;
    relrank0:=0; relrank1:=0; relrank2:=0;
    unique_marked:={}; unique_target_marked:={}; unique_halves:={};
    target_slice_points:={}; target_AB_pairs:={};
    order_cache:=AssociativeArray(); div_cache:=AssociativeArray();
    collision_signatures:=AssociativeArray();

    for A0 in K do
      if A0 eq 0 then continue; end if;
      for B0 in K do
        if B0 eq 0 then continue; end if;
        C0:=1/(A0*B0);
        R0:=A0^2+B0^2+C0^2-3;
        S0:=1/A0^2+1/B0^2+1/C0^2-3;
        if not IsSquare(R0) or not IsSquare(S0) then continue; end if;
        double_square_pairs+:=1;
        rz:=R0 eq 0; sz:=S0 eq 0;
        if rz then rho_zero_pairs+:=1; end if;
        if sz then sigma_zero_pairs+:=1; end if;
        if rz or sz then
            sheet_boundary_pairs+:=1;
            fprintf msk,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t-1\t0\t0\t0\t0\t0\n",
                    p,FInt(A0),FInt(B0),FInt(R0),FInt(S0),
                    rz select 1 else 0,sz select 1 else 0;
            continue;
        end if;

        r0:=SquareRoot(R0); s0:=SquareRoot(S0);
        // In the x=s^2 X normalization the four nonzero branch parameters
        // have squares A^2,B^2,C^2,(sigma/rho)^2=S/R.
        root2:=[A0^2,B0^2,C0^2,S0/R0];
        labels:=["A=B","A=C","A=D","B=C","B=D","C=D"];
        ij:=[<1,2>,<1,3>,<1,4>,<2,3>,<2,4>,<3,4>];
        equalities:=[labels[k]:k in [1..6]|root2[ij[k][1]] eq root2[ij[k][2]]];
        sig:=#equalities eq 0 select "smooth" else Join(equalities,"+");
        Bump(~collision_signatures,sig);
        rhos:=[r0,-r0]; sigmas:=[s0,-s0];
        base_collision:=false; base_smooth_sheets:=0;
        base_open_presentations:=0; base_target_presentations:=0;
        base_target_halves:=0; base_half_boundary:=0;

        for rho in rhos do for sigma in sigmas do
            affine_sheet_points+:=1;
            ss:=1/(2*rho); tt:=ss^2; uu:=2*tt;
            a0:=ss*A0; b0:=ss*B0; c0:=ss*C0;
            d0:=2*ss^2*sigma; vv:=3*tt^2+d0^2/2;
            P0,f0,q0,h0:=DirectData(K,a0,b0,c0,d0,uu,tt,vv);
            if h0^2-f0-q0^3 ne 0 then
                contact_failures+:=1;
                continue;
            end if;
            if Discriminant(f0) eq 0 then
                collision_sheet_points+:=1; base_collision:=true;
                continue;
            end if;
            smooth_sheet_points+:=1; base_smooth_sheets+:=1;

            triples:=[<a0,b0,c0>,<a0,c0,b0>,<b0,c0,a0>];
            for pairing in [1..3] do
                tr:=triples[pairing];
                aa:=tr[1]; bb:=tr[2]; cc:=tr[3]; dd:=d0;
                ok,f,q,h,g4,L4,U,V,mc,nc,det,ell:=
                    MarkedD12Aux(K,aa,bb,cc,dd,uu,tt,vv);
                if not ok then
                    chart_failures+:=1;
                    continue;
                end if;
                open_presentations+:=1; base_open_presentations+:=1;
                mkey:=MarkedKey(f,U,V);
                Include(~unique_marked,mkey);

                // Exact order check for every smooth/open marked presentation.
                orderD:=0; group_divisible:=false;
                if IsDefined(order_cache,mkey) then
                    orderD:=order_cache[mkey];
                    group_divisible:=div_cache[mkey];
                else
                    try
                        Jtmp:=Jacobian(HyperellipticCurve(f));
                        Dtmp:=Jtmp![U,V];
                        orderD:=Order(Dtmp);
                        if orderD eq 12 then
                            group_divisible:=IsDivisibleBy2Finite(Jtmp,Dtmp);
                        end if;
                        finite_group_decompositions+:=1;
                    catch e
                        orderD:=0; group_divisible:=false;
                    end try;
                    order_cache[mkey]:=orderD;
                    div_cache[mkey]:=group_divisible;
                end if;
                order_tests+:=1;
                if orderD eq 12 then order12_count+:=1;
                else order_other_count+:=1;
                end if;
                if orderD eq 12 then divisibility_tests+:=1; end if;

                okH,E1,E0,Spoly,ellH:=HalvingData(U,V,f);
                if not okH then continue; end if;
                d11:=Derivative(E1,1); d12:=Derivative(E1,2);
                d21:=Derivative(E0,1); d22:=Derivative(E0,2);
                class_has_half:=false; class_affine_halves:=0;
                P:=Parent(f); X:=P.1;
                for mm in K do for nn in K do
                    vals:=[mm,nn];
                    if Evaluate(E1,vals) ne 0 or Evaluate(E0,vals) ne 0 then
                        continue;
                    end if;
                    s4v:=Evaluate(Coefficient(Spoly,4),vals);
                    if s4v eq 0 then
                        half_boundary_points+:=1; base_half_boundary+:=1;
                        fprintf bdy,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\thalving_s4_zero\n",
                                p,FInt(A0),FInt(B0),FInt(C0),FInt(rho),
                                FInt(sigma),pairing,FInt(mm),FInt(nn);
                        continue;
                    end if;
                    s3v:=Evaluate(Coefficient(Spoly,3),vals);
                    s2v:=Evaluate(Coefficient(Spoly,2),vals);
                    G:=X^2+s3v/(2*s4v)*X+
                       (4*s4v*s2v-s3v^2)/(8*s4v^2);
                    ellv:=EvalCoeffPolynomial(ellH,vals,P);
                    W:=(-ellv) mod G;
                    if (W^2-f) mod G ne 0 then
                        error "square-quartic equations produced an invalid half";
                    end if;
                    ordH:=0; verified:=false;
                    if orderD eq 12 then
                        try
                            J:=Jacobian(HyperellipticCurve(f));
                            D:=J![U,V];
                            HH:=J![G,W];
                            verified:=2*HH eq D;
                            ordH:=Order(HH);
                        catch e
                            verified:=false;
                        end try;
                    end if;
                    if not verified or ordH ne 24 then
                        error "finite Jacobian half verification failed";
                    end if;
                    verified_halves+:=1;
                    class_has_half:=true;
                    class_affine_halves+:=1;
                    target_halves+:=1; base_target_halves+:=1;

                    rmat:=Rank(Matrix(K,2,2,[Evaluate(d11,vals),
                                              Evaluate(d12,vals),
                                              Evaluate(d21,vals),
                                              Evaluate(d22,vals)]));
                    if rmat eq 0 then relrank0+:=1;
                    elif rmat eq 1 then relrank1+:=1;
                    else relrank2+:=1;
                    end if;

                    Include(~unique_halves,HalfKey(f,U,V,G,W));
                    Include(~target_slice_points,
                            Sprint(<FInt(A0),FInt(B0),FInt(rho),FInt(sigma)>));
                    Include(~target_AB_pairs,Sprint(<FInt(A0),FInt(B0)>));
                    fprintf pts,
                      "%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                      p,FInt(A0),FInt(B0),FInt(C0),FInt(rho),FInt(sigma),
                      pairing,FInt(aa),FInt(bb),FInt(cc),FInt(dd),
                      FInt(mc),FInt(nc),FInt(Coefficient(U,1)),
                      FInt(Coefficient(U,0)),FInt(Coefficient(V,1)),
                      FInt(Coefficient(V,0)),FInt(mm),FInt(nn),
                      FInt(Coefficient(G,1)),FInt(Coefficient(G,0)),
                      FInt(Coefficient(W,1)),FInt(Coefficient(W,0)),
                      orderD,ordH,rmat;
                end for; end for;
                if class_has_half and not group_divisible then
                    error "affine square-quartic half disagrees with finite-group divisibility";
                end if;
                if group_divisible then
                    // The curve has full rational J[2], so a divisible class
                    // has exactly 16 group-theoretic halves.  Those absent
                    // from this affine (M,N) chart live on its projective
                    // boundary; retain their exact multiplicity separately.
                    missing:=16-class_affine_halves;
                    if missing lt 0 then
                        error "more than 16 affine halves for one marked class";
                    end if;
                    divisible_presentations+:=1;
                    target_presentations+:=1; base_target_presentations+:=1;
                    Include(~unique_target_marked,MarkedKey(f,U,V));
                    Include(~target_slice_points,
                            Sprint(<FInt(A0),FInt(B0),FInt(rho),FInt(sigma)>));
                    Include(~target_AB_pairs,Sprint(<FInt(A0),FInt(B0)>));
                    if missing gt 0 then
                        projective_MN_halves+:=missing;
                        fprintf bdy,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t-1\t-1\tprojective_MN_missing_%o\n",
                                p,FInt(A0),FInt(B0),FInt(C0),FInt(rho),
                                FInt(sigma),pairing,missing;
                    end if;
                end if;
            end for;
        end for; end for;

        fprintf msk,"%o\t%o\t%o\t%o\t%o\t0\t0\t%o\t%o\t%o\t%o\t%o\t%o\n",
                p,FInt(A0),FInt(B0),FInt(R0),FInt(S0),
                base_collision select 1 else 0,base_smooth_sheets,
                base_open_presentations,base_target_presentations,
                base_target_halves,base_half_boundary;
      end for;
    end for;

    delete pts; delete bdy; delete msk;
    print "FINITE_EXHAUSTIVE",p,
          "scope","A*B*rho*sigma*disc*CRTdet*s4 != 0",
          "torus_AB_pairs",torus_pairs,
          "double_square_AB_pairs",double_square_pairs,
          "rho_zero_pairs",rho_zero_pairs,
          "sigma_zero_pairs",sigma_zero_pairs,
          "sheet_boundary_AB_pairs",sheet_boundary_pairs,
          "affine_slice_sheet_points",affine_sheet_points,
          "collision_sheet_points",collision_sheet_points,
          "smooth_slice_sheet_points",smooth_sheet_points,
          "open_marked_presentations",open_presentations,
          "CRT_chart_failures",chart_failures,
          "contact_failures",contact_failures;
    print "FINITE_ORDERS",p,
          "tests",order_tests,"order12",order12_count,"other",order_other_count,
          "divisibility_tests",divisibility_tests,
          "finite_group_decompositions",finite_group_decompositions,
          "divisible_presentations",divisible_presentations,
          "unique_marked_classes",#unique_marked;
    print "FINITE_TARGETS",p,
          "target_AB_pairs",#target_AB_pairs,
          "target_slice_points",#target_slice_points,
          "target_presentations",target_presentations,
          "unique_target_marked",#unique_target_marked,
          "target_halves",target_halves,
          "unique_halves",#unique_halves,
          "verified_halves",verified_halves,
          "relative_rank_counts",[relrank0,relrank1,relrank2],
          "finite_s4_boundary_points",half_boundary_points,
          "projective_MN_halves",projective_MN_halves,
          "all_group_halves",target_halves+projective_MN_halves;
    print "FINITE_FILES",p,points_name,boundary_name,mask_name;
    print "FINITE_PROJECTIVE_SCOPE",p,
          "A=0, B=0, A=infinity, B=infinity, rho=0, sigma=0, and s4=0 are boundary strata, not affine obstructions";
    print "FINITE_BRANCH_SIGNATURES",p,
          Sort([<k,collision_signatures[k]>:k in Keys(collision_signatures)]);
end procedure;

print "TARGET_22224_QSQUARE_CRT_SIEVE_START";
print "CONFIG","PrimeList",PrimeList,"output_stem",output_stem;
for p in PrimeList do
    EnumeratePrime(p,output_stem);
end for;
print "TARGET_22224_QSQUARE_CRT_SIEVE_DONE";
UnsetLogFile();
quit;
