//////////////////////////////////////////////////////////////////////
// Family-level Richelot sieve for target [2,2,4,8].
//
// Stage A (reverse generation): enumerate primitive square-branch curves
//
//   y^2=x(x+a^2)(x+b^2)(x+c^2)(x+d^2),
//
// enumerate all 15 rational Richelot neighbors, and retain the distinct
// one-split neighbors with exact torsion [2,4,8] and a strict root-power
// geometric-simplicity certificate.  This efficiently recovers a bank of
// sources without a quadratic-height box on the M_1(8,4) cover.
//
// Stage B (forward target sieve): enumerate every rational Richelot neighbor
// of every retained source.  A [2,2,4,8] target must have six rational
// Weierstrass points and 2^7 dividing every good-reduction Jacobian order.
// These two cheap filters precede exact TorsionSubgroup calls.
//
// Typical run from torsion_jac:
//   magma -b tuple_bound:=10 max_bases:=0 max_sources:=50 \
//     write_log:=true code/target_2248_family_richelot_sieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned tuple_bound then tuple_bound:=10;
elif Type(tuple_bound) eq MonStgElt then tuple_bound:=StringToInteger(tuple_bound); end if;
if not assigned max_bases then max_bases:=0;
elif Type(max_bases) eq MonStgElt then max_bases:=StringToInteger(max_bases); end if;
if not assigned max_sources then max_sources:=50;
elif Type(max_sources) eq MonStgElt then max_sources:=StringToInteger(max_sources); end if;
if not assigned reduction_prime_count then reduction_prime_count:=2;
elif Type(reduction_prime_count) eq MonStgElt then reduction_prime_count:=StringToInteger(reduction_prime_count); end if;
if not assigned exact_all_full then exact_all_full:=false;
elif Type(exact_all_full) eq MonStgElt then exact_all_full:=exact_all_full in {"true","True","1","yes"}; end if;
if not assigned include_controls then include_controls:=true;
elif Type(include_controls) eq MonStgElt then include_controls:=include_controls in {"true","True","1","yes"}; end if;
if not assigned write_log then write_log:=true;
elif Type(write_log) eq MonStgElt then write_log:=write_log in {"true","True","1","yes"}; end if;
if not assigned log_file then log_file:="results/target_2248_family_richelot_sieve.log"; end if;
if write_log then SetLogFile(log_file : Overwrite:=true); end if;

Q:=Rationals(); Z:=Integers(); P<x>:=PolynomialRing(Q);
prime_list:=[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
             67,71,73,79,83,89,97,101,103,107,109,113,127,131];

function BoolString(b) return b select "true" else "false"; end function;
function InvOrder(inv) return #inv eq 0 select 1 else &*inv; end function;

function IntegralSquareScale(f)
    d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return P!(d^2*f);
end function;

function NormalizeCurve(C)
    f,h:=HyperellipticPolynomials(C);
    F:=h eq 0 select P!f else P!(h^2+4*f);
    D:=SimplifiedModel(HyperellipticCurve(IntegralSquareScale(F)));
    fD,hD:=HyperellipticPolynomials(D);
    return D,P!fD,P!hD;
end function;

function CompletedPolynomial(C)
    f,h:=HyperellipticPolynomials(C);
    return h eq 0 select P!f else P!(h^2+4*f);
end function;

function AbsoluteInvariants(C)
    return [Q!a:a in G2Invariants(C)];
end function;

function SameQCurve(C,D)
    try ok,mp:=IsIsomorphic(C,D); return ok;
    catch e return false; end try;
end function;

function ExistingCurve(C,ai,curves,ais)
    for i in [1..#curves] do
        if ais[i] eq ai and SameQCurve(C,curves[i]) then return i; end if;
    end for;
    return 0;
end function;

function FactorDegrees(C)
    F:=CompletedPolynomial(C);
    return Sort(&cat[[Degree(z[1]):j in [1..z[2]]]:z in Factorization(F)]);
end function;

function OneSplitPattern(C)
    F:=CompletedPolynomial(C); ds:=FactorDegrees(C);
    if Degree(F) eq 5 then return ds eq [1,1,1,2]; end if;
    return Degree(F) eq 6 and ds eq [1,1,1,1,2];
end function;

function FullWeierstrassPattern(C)
    F:=CompletedPolynomial(C); ds:=FactorDegrees(C);
    if Degree(F) eq 5 then return ds eq [1,1,1,1,1]; end if;
    return Degree(F) eq 6 and ds eq [1,1,1,1,1,1];
end function;

function TwoAdicReductionGate(C,need_v2)
    F:=CompletedPolynomial(C); used:=0; g:=0; rows:=[];
    for pp in prime_list do
        try
            Fp:=ChangeRing(F,GF(pp));
            if Degree(Fp) ne Degree(F) or Discriminant(Fp) eq 0 then continue; end if;
            Lp:=LPolynomial(ChangeRing(C,GF(pp)));
            n:=Z!Evaluate(Lp,1);
            g:=used eq 0 select n else GCD(g,n);
            used+:=1; Append(~rows,<pp,n,Valuation(n,2)>);
            if Valuation(g,2) lt need_v2 then return false,g,rows; end if;
            if used ge reduction_prime_count then return true,g,rows; end if;
        catch e
            continue;
        end try;
    end for;
    return false,g,rows;
end function;

function StrictRootPowerCertificate(C)
    F:=CompletedPolynomial(C);
    for pp in prime_list do
        try
            Fp:=ChangeRing(F,GF(pp));
            if Degree(Fp) ne Degree(F) or Discriminant(Fp) eq 0 then continue; end if;
            Lp:=LPolynomial(ChangeRing(C,GF(pp)));
            cs:=[Z!Coefficient(Lp,i):i in [0..4]];
            S<T>:=PolynomialRing(Q);
            chi:=&+[Q!cs[5-i]*T^i:i in [0..4]];
            if not IsIrreducible(chi) then continue; end if;
            K<pi>:=NumberField(chi);
            ds:=[Degree(MinimalPolynomial(pi^n)):n in [2..12]];
            if &and[d eq 4:d in ds] then return true,pp,chi,ds; end if;
        catch e
            continue;
        end try;
    end for;
    S<T>:=PolynomialRing(Q);
    return false,0,S!0,[];
end function;

function Bracket(A,B) return Derivative(A)*B-A*Derivative(B); end function;

function CoefficientMatrixDeterminant(A,B,C)
    M:=Matrix(Q,3,3,[Coefficient(A,i):i in [0..2]] cat
                     [Coefficient(B,i):i in [0..2]] cat
                     [Coefficient(C,i):i in [0..2]]);
    return Determinant(M);
end function;

function Pairings4()
    return [[<1,2>,<3,4>],[<1,3>,<2,4>],[<1,4>,<2,3>]];
end function;

// For a target edge, identify the built-in codomain with its explicit
// Galois-stable quadratic splitting.  One factor may be linear when it pairs
// a finite root with infinity in an odd model.
procedure PrintKernelIdentification(Csource,Ctarget,builtin_index)
    F:=CompletedPolynomial(Csource); fac:=Factorization(F);
    lins:=[z[1]:z in fac|Degree(z[1]) eq 1];
    quads:=[z[1]:z in fac|Degree(z[1]) eq 2];
    if #quads ne 1 then
        print "TARGET_KERNEL_IDENTIFICATION_UNAVAILABLE","factorization",fac;
        return;
    end if;
    triples:=[]; labels:=[];
    if Degree(F) eq 5 and #lins eq 3 then
        for k in [1..3] do
            ij:=[i:i in [1..3]|i ne k];
            Append(~triples,<quads[1],lins[ij[1]]*lins[ij[2]],lins[k]>);
            Append(~labels,Sprintf("Q | L%o*L%o | L%o*infinity",ij[1],ij[2],k));
        end for;
    elif Degree(F) eq 6 and #lins eq 4 then
        for pairing in Pairings4() do
            Append(~triples,<quads[1],
                    lins[pairing[1][1]]*lins[pairing[1][2]],
                    lins[pairing[2][1]]*lins[pairing[2][2]]>);
            Append(~labels,Sprintf("Q | L%o*L%o | L%o*L%o",
                pairing[1][1],pairing[1][2],pairing[2][1],pairing[2][2]));
        end for;
    else
        print "TARGET_KERNEL_IDENTIFICATION_UNAVAILABLE","factorization",fac;
        return;
    end if;
    aiT:=AbsoluteInvariants(Ctarget);
    for k in [1..#triples] do
        A:=triples[k][1]; B:=triples[k][2]; D:=triples[k][3];
        delta:=CoefficientMatrixDeterminant(A,B,D);
        g0:=delta*Bracket(B,D)*Bracket(D,A)*Bracket(A,B);
        if delta eq 0 or Discriminant(g0) eq 0 then continue; end if;
        Cp,fp,hp:=NormalizeCurve(HyperellipticCurve(g0));
        Cm,fm,hm:=NormalizeCurve(HyperellipticCurve(-g0));
        plus:=AbsoluteInvariants(Cp) eq aiT and SameQCurve(Cp,Ctarget);
        minus:=AbsoluteInvariants(Cm) eq aiT and SameQCurve(Cm,Ctarget);
        if plus or minus then
            print "TARGET_KERNEL_IDENTIFICATION","builtin_index",builtin_index,
                  "kernel_index",k,"label",labels[k],"factors",[A,B,D],
                  "delta",delta,"formula_sign",plus select 1 else -1,
                  "manual_g0",g0;
        end if;
    end for;
end procedure;

print "TARGET_2248_FAMILY_RICHELOT_SIEVE";
print "CONFIG","tuple_bound",tuple_bound,"max_bases",max_bases,
      "max_sources",max_sources,"reduction_prime_count",reduction_prime_count,
      "exact_all_full",BoolString(exact_all_full),
      "include_controls",BoolString(include_controls);

base_seen_curves:=[]; base_seen_ai:=[];
candidate_seen_curves:=[]; candidate_seen_ai:=[];
source_curves:=[]; source_fs:=[]; source_hs:=[]; source_ai:=[];
source_origin_tuple:=[]; source_origin_edge:=[];
source_cert_p:=[]; source_cert_chi:=[]; source_cert_degrees:=[];

// Exact unique polynomials represented by the eight sign/reciprocal hits in
// data/m18_m14_one_split_h20.txt.  They are positive controls for source
// recovery and are Q-isomorphism-deduplicated below.
if include_controls then
    control_polynomials:=[
      1569052980450*x^5 + 17480455510401*x^4 + 53939339809344*x^3
        + 65460742603776*x^2 + 27862813900800*x,
      7061463847622250*x^5 + 104632219276049025*x^4
        + 135735215960638800*x^3 + 188573481843278400*x^2
        + 51200550567936000*x,
      484884627849216000*x^5 + 1785847643268710400*x^4
        + 1285453358353612800*x^3 + 990898616165990400*x^2
        + 66874188496896000*x,
      879609302220800*x^5 + 2066549284274176*x^4
        + 1702826757586944*x^3 + 551845600690176*x^2
        + 49533891379200*x
    ];
    control_curves:=[]; control_ai:=[];
    for k in [1..#control_polynomials] do
        D,fD,hD:=NormalizeCurve(HyperellipticCurve(control_polynomials[k]));
        aiD:=AbsoluteInvariants(D);
        old:=ExistingCurve(D,aiD,control_curves,control_ai);
        if old ne 0 then
            print "CONTROL_DUPLICATE",k,"existing_control",old;
            continue;
        end if;
        Append(~control_curves,D); Append(~control_ai,aiD);
        G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
        simple,pcert,chi,degrees:=StrictRootPowerCertificate(D);
        print "CONTROL_SOURCE",k,"unique_control",#control_curves,
              "factor_degrees",FactorDegrees(D),"torsion",inv,"order",#G,
              "root_power",BoolString(simple),"p",pcert,"chi",chi,
              "degrees",degrees,"f",fD,"h",hD,"absolute_G2",aiD;
        if inv eq [2,4,8] and simple and
           ExistingCurve(D,aiD,source_curves,source_ai) eq 0 then
            Append(~source_curves,D); Append(~source_fs,fD); Append(~source_hs,hD);
            Append(~source_ai,aiD); Append(~source_origin_tuple,[0,0,0,k]);
            Append(~source_origin_edge,0); Append(~source_cert_p,pcert);
            Append(~source_cert_chi,chi); Append(~source_cert_degrees,degrees);
            print "CONTROL_SOURCE_ACCEPT",#source_curves,"control",k;
        end if;
    end for;
    print "CONTROL_SUMMARY","input_polynomials",#control_polynomials,
          "unique_Q_curves",#control_curves,"accepted_simple_248",#source_curves;
end if;

tuples_checked:=0; primitive_tuples:=0; bases_used:=0; base_duplicates:=0;
reverse_edges:=0; reverse_jacobians:=0; reverse_one_split:=0;
reverse_reduction_pass:=0; reverse_candidate_duplicates:=0;
reverse_exact_tests:=0; reverse_exact_248:=0; reverse_simple_248:=0;

stop_generation:=max_sources gt 0 and #source_curves ge max_sources;
for a in [1..tuple_bound-3] do
  if stop_generation then break; end if;
  for b in [a+1..tuple_bound-2] do
    if stop_generation then break; end if;
    for c in [b+1..tuple_bound-1] do
      if stop_generation then break; end if;
      for d in [c+1..tuple_bound] do
        tuples_checked+:=1;
        if GCD(GCD(a,b),GCD(c,d)) ne 1 then continue; end if;
        primitive_tuples+:=1;
        f0:=x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
        C0,f0n,h0n:=NormalizeCurve(HyperellipticCurve(f0));
        ai0:=AbsoluteInvariants(C0);
        if ExistingCurve(C0,ai0,base_seen_curves,base_seen_ai) ne 0 then
            base_duplicates+:=1; continue;
        end if;
        if max_bases gt 0 and bases_used ge max_bases then
            stop_generation:=true; break;
        end if;
        Append(~base_seen_curves,C0); Append(~base_seen_ai,ai0);
        bases_used+:=1;
        if bases_used mod 10 eq 0 then
            print "GENERATION_PROGRESS","bases",bases_used,"tuples",tuples_checked,
                  "edges",reverse_edges,"one_split",reverse_one_split,
                  "exact_tests",reverse_exact_tests,"sources",#source_curves;
        end if;
        try
            Rs:=RichelotIsogenousSurfaces(Jacobian(C0));
            print "BASE",bases_used,"tuple",[a,b,c,d],"richelot_count",#Rs;
            for j in [1..#Rs] do
                reverse_edges+:=1;
                if Type(Rs[j]) ne JacHyp then continue; end if;
                reverse_jacobians+:=1;
                D,fD,hD:=NormalizeCurve(Curve(Rs[j]));
                if not OneSplitPattern(D) then continue; end if;
                reverse_one_split+:=1;
                gate,g,redrows:=TwoAdicReductionGate(D,6);
                if not gate then continue; end if;
                reverse_reduction_pass+:=1;
                aiD:=AbsoluteInvariants(D);
                if ExistingCurve(D,aiD,candidate_seen_curves,candidate_seen_ai) ne 0 then
                    reverse_candidate_duplicates+:=1; continue;
                end if;
                Append(~candidate_seen_curves,D); Append(~candidate_seen_ai,aiD);
                reverse_exact_tests+:=1;
                G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
                print "SOURCE_CANDIDATE","base",bases_used,"tuple",[a,b,c,d],
                      "edge",j,"factor_degrees",FactorDegrees(D),
                      "reduction_gcd",g,"reduction_rows",redrows,
                      "torsion",inv,"order",#G;
                if inv ne [2,4,8] then continue; end if;
                reverse_exact_248+:=1;
                simple,pcert,chi,degrees:=StrictRootPowerCertificate(D);
                print "SOURCE_248_CERT","simple",BoolString(simple),"p",pcert,
                      "chi",chi,"degrees",degrees;
                if not simple then continue; end if;
                reverse_simple_248+:=1;
                if ExistingCurve(D,aiD,source_curves,source_ai) ne 0 then continue; end if;
                Append(~source_curves,D); Append(~source_fs,fD); Append(~source_hs,hD);
                Append(~source_ai,aiD); Append(~source_origin_tuple,[a,b,c,d]);
                Append(~source_origin_edge,j); Append(~source_cert_p,pcert);
                Append(~source_cert_chi,chi); Append(~source_cert_degrees,degrees);
                si:=#source_curves;
                print "SOURCE_ACCEPT",si,"origin_tuple",[a,b,c,d],"origin_edge",j,
                      "f",fD,"h",hD,"absolute_G2",aiD,
                      "root_power_p",pcert,"root_power_chi",chi;
                if max_sources gt 0 and #source_curves ge max_sources then
                    stop_generation:=true; break;
                end if;
            end for;
        catch e
            print "BASE_ERROR",bases_used,"tuple",[a,b,c,d],e`Object;
        end try;
        if stop_generation then break; end if;
      end for;
    end for;
  end for;
end for;

print "GENERATION_SUMMARY","tuples_checked",tuples_checked,
      "primitive_tuples",primitive_tuples,"bases_used",bases_used,
      "base_duplicates",base_duplicates,"edges",reverse_edges,
      "jacobians",reverse_jacobians,"one_split",reverse_one_split,
      "reduction_pass",reverse_reduction_pass,
      "candidate_duplicates",reverse_candidate_duplicates,
      "exact_tests",reverse_exact_tests,"exact_248",reverse_exact_248,
      "simple_248",reverse_simple_248,"accepted_sources",#source_curves;

forward_edges:=0; forward_jacobians:=0; forward_full:=0;
forward_reduction_pass:=0; forward_exact_tests:=0; target_hits:=0;
neighbor_seen_curves:=[]; neighbor_seen_ai:=[]; neighbor_duplicates:=0;

for si in [1..#source_curves] do
    C:=source_curves[si];
    print "SOURCE_FORWARD",si,"origin_tuple",source_origin_tuple[si],
          "origin_edge",source_origin_edge[si],"f",source_fs[si],"h",source_hs[si],
          "cert_p",source_cert_p[si],"cert_chi",source_cert_chi[si],
          "cert_degrees",source_cert_degrees[si];
    try
        Rs:=RichelotIsogenousSurfaces(Jacobian(C));
        print "SOURCE_FORWARD_RICHELOT_COUNT",si,#Rs;
        for j in [1..#Rs] do
            forward_edges+:=1;
            if Type(Rs[j]) ne JacHyp then
                print "FORWARD_EDGE",si,j,"type",Type(Rs[j]); continue;
            end if;
            forward_jacobians+:=1;
            D,fD,hD:=NormalizeCurve(Curve(Rs[j]));
            ds:=FactorDegrees(D); full:=FullWeierstrassPattern(D);
            aiD:=AbsoluteInvariants(D);
            old:=ExistingCurve(D,aiD,neighbor_seen_curves,neighbor_seen_ai);
            if old ne 0 then neighbor_duplicates+:=1;
            else Append(~neighbor_seen_curves,D); Append(~neighbor_seen_ai,aiD); end if;
            print "FORWARD_EDGE",si,j,"type JacHyp","factor_degrees",ds,
                  "full_weierstrass",BoolString(full),"global_neighbor_index",old;
            if not full then continue; end if;
            forward_full+:=1;
            gate,g,redrows:=TwoAdicReductionGate(D,7);
            print "FORWARD_FULL_GATE",si,j,"pass",BoolString(gate),
                  "gcd",g,"rows",redrows;
            if not gate and not exact_all_full then continue; end if;
            if gate then forward_reduction_pass+:=1; end if;
            forward_exact_tests+:=1;
            G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
            print "FORWARD_EXACT",si,j,"torsion",inv,"order",#G,
                  "f",fD,"h",hD,"absolute_G2",aiD;
            if inv eq [2,2,4,8] then
                target_hits+:=1;
                simpleT,pT,chiT,dsT:=StrictRootPowerCertificate(D);
                print "TARGET_2248_HIT","source",si,"edge",j,
                      "torsion",inv,"source_root_power_p",source_cert_p[si],
                      "source_root_power_chi",source_cert_chi[si],
                      "target_root_power",BoolString(simpleT),"target_p",pT,
                      "target_chi",chiT,"target_degrees",dsT,
                      "target_f",fD,"target_h",hD;
                PrintKernelIdentification(C,D,j);
            end if;
        end for;
    catch e
        print "SOURCE_FORWARD_ERROR",si,e`Object;
    end try;
end for;

print "FORWARD_SUMMARY","sources",#source_curves,"edges",forward_edges,
      "jacobians",forward_jacobians,"full_weierstrass",forward_full,
      "reduction_pass",forward_reduction_pass,"exact_tests",forward_exact_tests,
      "neighbor_duplicates",neighbor_duplicates,
      "unique_neighbors",#neighbor_seen_curves,"target_hits",target_hits;
print "TARGET_2248_FAMILY_RICHELOT_SIEVE_DONE";
quit;
