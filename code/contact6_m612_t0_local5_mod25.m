//////////////////////////////////////////////////////////////////////
// Exact first blow-up / mod-25 diagnostic for the integral T0 cover.
//
// Every mod-5 point x0 of the four-equation pullback is lifted as
// x=x0+5*x1.  The mod-25 equations are the exact affine linear system
//
//   J(x0)*x1 = -F(x0)/5  (mod 5).
//
// We then require every open factor which vanished at x0 to have exact
// valuation one at x.  Such a lift is a genuine first-order cone away from
// all base-discriminant and tangent/contact divisors.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers(); k:=GF(5);
R<s,omega,m,L,U,nu>:=PolynomialRing(Z,6,"grevlex");

a:=s*omega^2-3; b:=s-3; W:=m^2; M:=L^2;
c5:=2*s; c4:=b^2+2*a-15; c3:=2*a*b+22;
c2:=a^2+2*b-15; c1:=2*(a+3);
N:=c2-omega*c4+omega*W;
R8:=8*s*(c3-4*s*omega)-(c4-W)^2;
Fhalf:=R8^2-256*s^2*W*N;
B3:=c5*M+3*U;
Delta3:=4*c4*M+12*(U^2+nu^2)-B3^2;
F3:=B3*Delta3+16*nu^3-8*c3*M-8*U^3-48*U*nu^2;
F2:=Delta3^2+64*B3*nu^3-64*c2*M
    -192*(U^2*nu^2+nu^4);
F1:=Delta3*nu^3-4*c1*M-12*U*nu^4;
eqs:=[Fhalf,F1,F2,F3];
derivs:=[[Derivative(f,j):j in [1..6]]:f in eqs];

DB:=(a-3)^2-8*(b+3);
DC:=(b-3)^2-8*(a+3);
base_names:=["b3","a3","RR","DB","DC"];
base_forms:=[s,s*omega^2,a+b+2,DB,DC];
deep_names:=["m","L","nu","Up","Um","omega_pm_1"];
deep_forms:=[m,L,nu,U-2*nu,U+2*nu,omega^2-1];
all_names:=base_names cat deep_names;
all_forms:=base_forms cat deep_forms;

function Mod5Z(z)
    return k!(z mod 5);
end function;

function VanishingSet(forms,pt)
    return [i:i in [1..#forms]|(Z!Evaluate(forms[i],pt)) mod 5 eq 0];
end function;

residue_points:=0; mod25_lifts:=0; escape_base:=0;
escape_deep:=0; escape_all:=0;
by_cone:=AssociativeArray(); by_rank:=AssociativeArray();
by_cone_rank:=AssociativeArray();
samples:=[]; smooth_samples:=[]; smooth_goodbase_samples:=[];
forced_rank4_pairs:={@ @}; forced_rank4_full:={@ @};
forced_pair_counts:=AssociativeArray();
forced_rank3_pairs:={@ @}; forced_rank3_full:={@ @};
forced_rank3_pair_counts:=AssociativeArray(); forced_rank3_samples:=[];

for s0 in [0..4] do for o0 in [0..4] do for m0 in [0..4] do
for l0 in [0..4] do for u0 in [0..4] do for n0 in [0..4] do
    x0:=[s0,o0,m0,l0,u0,n0];
    f0:=[Z!Evaluate(f,x0):f in eqs];
    if not &and[z mod 5 eq 0:z in f0] then continue; end if;
    residue_points+:=1;
    A:=Matrix(k,4,6,
        &cat[[Mod5Z(Evaluate(derivs[i][j],x0)):j in [1..6]]
              :i in [1..4]]);
    rhs:=Vector(k,[Mod5Z(-(f0[i] div 5)):i in [1..4]]);
    rk:=Rank(A);
    bzero:=VanishingSet(base_forms,x0);
    dzero:=VanishingSet(deep_forms,x0);
    for s1 in [0..4] do for o1 in [0..4] do for m1 in [0..4] do
    for l1 in [0..4] do for u1 in [0..4] do for n1 in [0..4] do
        x1:=[s1,o1,m1,l1,u1,n1];
        if Vector(k,x1)*Transpose(A) ne rhs then continue; end if;
        mod25_lifts+:=1;
        x:=[x0[i]+5*x1[i]:i in [1..6]];
        // Since nonvanishing residues are already units, it is enough to
        // test exact escape for the factors which vanish at x0.
        eb:=&and[(Z!Evaluate(base_forms[i],x)) mod 25 ne 0:i in bzero];
        ed:=&and[(Z!Evaluate(deep_forms[i],x)) mod 25 ne 0:i in dzero];
        if eb then escape_base+:=1; end if;
        if ed then escape_deep+:=1; end if;
        if not (eb and ed) then continue; end if;
        escape_all+:=1;
        bname:=#bzero eq 0 select "goodbase"
               else Join([base_names[i]:i in bzero],"+");
        dname:=#dzero eq 0 select "deepopen"
               else Join([deep_names[i]:i in dzero],"+");
        key:=bname cat "__" cat dname;
        if IsDefined(by_cone,key) then by_cone[key]+:=1;
        else by_cone[key]:=1; end if;
        if IsDefined(by_rank,rk) then by_rank[rk]+:=1;
        else by_rank[rk]:=1; end if;
        keyrk:=key cat "__rank" cat IntegerToString(rk);
        if IsDefined(by_cone_rank,keyrk) then by_cone_rank[keyrk]+:=1;
        else by_cone_rank[keyrk]:=1; end if;
        if #samples lt 4 then
            vals:=[(Z!Evaluate(all_forms[i],x)) mod 25:i in [1..#all_forms]];
            Append(~samples,Sprint(<x0,x1,x,key,rk,vals>));
        end if;
        if rk eq 4 and #smooth_samples lt 4 then
            vals:=[(Z!Evaluate(all_forms[i],x)) mod 25:i in [1..#all_forms]];
            Append(~smooth_samples,Sprint(<x0,x1,x,key,vals>));
        end if;
        if rk eq 4 and #bzero eq 0 and #smooth_goodbase_samples lt 4 then
            vals:=[(Z!Evaluate(all_forms[i],x)) mod 25:i in [1..#all_forms]];
            Append(~smooth_goodbase_samples,Sprint(<x0,x1,x,key,vals>));
        end if;
        if rk eq 4 and #bzero gt 0 then
            pair:=<x[2] mod 25,x[6] mod 25>;
            Include(~forced_rank4_pairs,pair);
            Include(~forced_rank4_full,<x[1] mod 25,x[2] mod 25,
                                      x[3] mod 25,x[4] mod 25,
                                      x[5] mod 25,x[6] mod 25>);
            pairkey:=Sprint(pair);
            if IsDefined(forced_pair_counts,pairkey) then
                forced_pair_counts[pairkey]+:=1;
            else
                forced_pair_counts[pairkey]:=1;
            end if;
        end if;
        if rk eq 3 and #bzero gt 0 then
            pair:=<x[2] mod 25,x[6] mod 25>;
            Include(~forced_rank3_pairs,pair);
            Include(~forced_rank3_full,<x[1] mod 25,x[2] mod 25,
                                      x[3] mod 25,x[4] mod 25,
                                      x[5] mod 25,x[6] mod 25>);
            pairkey:=Sprint(pair);
            if IsDefined(forced_rank3_pair_counts,pairkey) then
                forced_rank3_pair_counts[pairkey]+:=1;
            else
                forced_rank3_pair_counts[pairkey]:=1;
            end if;
            if #forced_rank3_samples lt 8 then
                exact_residuals:=[(Z!Evaluate(f,x)) mod 25:f in eqs];
                Append(~forced_rank3_samples,Sprint(<x0,x1,x,key,
                                                     exact_residuals>));
            end if;
        end if;
    end for; end for; end for; end for; end for; end for;
end for; end for; end for; end for; end for; end for;

print "CONTACT6_M612_T0_LOCAL5_MOD25";
print "RESIDUE_POINTS",residue_points,"MOD25_LIFTS",mod25_lifts;
print "ESCAPE_BASE",escape_base,"ESCAPE_DEEP",escape_deep,
      "ESCAPE_ALL",escape_all;
print "ESCAPE_ALL_BY_RESIDUE_RANK",
      Sort([<r,by_rank[r]>:r in Keys(by_rank)]);
print "ESCAPE_ALL_CONES",Sort([<x,by_cone[x]>:x in Keys(by_cone)]);
print "ESCAPE_ALL_CONE_RANKS",
      Sort([<x,by_cone_rank[x]>:x in Keys(by_cone_rank)]);
print "FACTOR_ORDER",all_names;
print "SAMPLES_x0_x1_xmod25_cone_rank_factor_values_mod25",samples;
print "SMOOTH_SAMPLES_x0_x1_xmod25_cone_factor_values_mod25",smooth_samples;
print "SMOOTH_GOODBASE_SAMPLES_x0_x1_xmod25_cone_factor_values_mod25",
      smooth_goodbase_samples;
print "FORCED_BOUNDARY_RANK4_OMEGA_NU_PAIRS_MOD25",
      Sort(Setseq(forced_rank4_pairs));
print "FORCED_BOUNDARY_RANK4_PAIR_COUNTS",
      Sort([<x,forced_pair_counts[x]>:x in Keys(forced_pair_counts)]);
print "FORCED_BOUNDARY_RANK4_FULL_POINTS_s_omega_m_L_U_nu_MOD25",
      Sort(Setseq(forced_rank4_full));
print "FORCED_BOUNDARY_RANK3_OMEGA_NU_PAIRS_MOD25",
      Sort(Setseq(forced_rank3_pairs));
print "FORCED_BOUNDARY_RANK3_PAIR_COUNTS",
      Sort([<x,forced_rank3_pair_counts[x]>:
            x in Keys(forced_rank3_pair_counts)]);
print "FORCED_BOUNDARY_RANK3_FULL_POINTS_s_omega_m_L_U_nu_MOD25",
      Sort(Setseq(forced_rank3_full));
print "FORCED_BOUNDARY_RANK3_EXACT_SAMPLES_x0_x1_x_cone_residuals",
      forced_rank3_samples;
print "CONTACT6_M612_T0_LOCAL5_MOD25_DONE";
quit;
