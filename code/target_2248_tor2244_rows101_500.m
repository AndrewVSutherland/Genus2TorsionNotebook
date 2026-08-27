//////////////////////////////////////////////////////////////////////
// Target [2,2,4,8]: rows 101--500 of the A(2,2,4,4) source bank.
//
// Each bank row defines
//     y^2 = x prod_{u in {a,b,c,d}} (x+u^2).
// We enumerate all rational Richelot neighbours, retain one-split
// quotients with exact torsion [2,4,8] and deduplicate them over Q.
// Every exact source is retained and classified as certified-simple,
// proven-product, or unresolved: failure of the strict root-power test
// alone is never treated as nonsimplicity.  Every full-Weierstrass
// rational Richelot neighbour is then subjected to an exact torsion call.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned first_row then first_row:=101;
elif Type(first_row) eq MonStgElt then first_row:=StringToInteger(first_row); end if;
if not assigned last_row then last_row:=500;
elif Type(last_row) eq MonStgElt then last_row:=StringToInteger(last_row); end if;
if not assigned bank_file then bank_file:="data/tor2244_bank.txt"; end if;
if not assigned log_file then
    log_file:="results/target_2248_tor2244_rows101_500.log";
end if;
if not assigned reduction_prime_count then reduction_prime_count:=2;
elif Type(reduction_prime_count) eq MonStgElt then
    reduction_prime_count:=StringToInteger(reduction_prime_count);
end if;
if not assigned write_log then write_log:=true;
elif Type(write_log) eq MonStgElt then
    write_log:=write_log in {"true","True","1","yes"};
end if;
if write_log then SetLogFile(log_file : Overwrite:=true); end if;

Q:=Rationals(); Z:=Integers(); P<x>:=PolynomialRing(Q);
prime_list:=[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
             67,71,73,79,83,89,97,101,103,107,109,113,127,131];

function BoolString(b) return b select "true" else "false"; end function;

function ReadTupleFile(filename)
    raw_rows:=Split(Read(filename),"\n"); out:=[];
    for raw in raw_rows do
        n:=#raw;
        if n gt 0 and raw[n] eq "\r" then n-:=1; end if;
        if n ge 2 and raw[1] eq "[" and raw[n] eq "]" then
            ss:=Split(raw[2..n-1],",");
            if #ss eq 4 then Append(~out,[StringToInteger(s):s in ss]); end if;
        end if;
    end for;
    return out;
end function;

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

// Necessary good-reduction divisibility gate: order 2^need_v2 must divide
// every good-reduction Jacobian order.
function TwoAdicReductionGate(C,need_v2)
    F:=CompletedPolynomial(C); used:=0; g:=0; redrows:=[];
    for pp in prime_list do
        try
            Fp:=ChangeRing(F,GF(pp));
            if Degree(Fp) ne Degree(F) or Discriminant(Fp) eq 0 then continue; end if;
            Lp:=LPolynomial(ChangeRing(C,GF(pp)));
            n:=Z!Evaluate(Lp,1);
            g:=used eq 0 select n else GCD(g,n);
            used+:=1; Append(~redrows,<pp,n,Valuation(n,2)>);
            if Valuation(g,2) lt need_v2 then return false,g,redrows; end if;
            if used ge reduction_prime_count then return true,g,redrows; end if;
        catch e
            continue;
        end try;
    end for;
    return false,g,redrows;
end function;

// An irreducible Frobenius polynomial whose roots retain degree four after
// powers 2,...,12 is the strict certificate used throughout this project.
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

// A rational Richelot quotient which is a Cartesian product proves that
// the source Jacobian is geometrically nonsimple.  Absence of such a
// rational product quotient proves nothing, so return only a positive
// certificate here.
function RationalProductCertificate(C)
    try
        Rs:=RichelotIsogenousSurfaces(Jacobian(C));
        for j in [1..#Rs] do
            if Sprint(Type(Rs[j])) eq "SetCart" then return true,j; end if;
        end for;
    catch e
        dummy:=0;
    end try;
    return false,0;
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

// If a target appears, identify its rational Richelot kernel explicitly.
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

rows:=ReadTupleFile(bank_file);
first:=Max(1,first_row); last:=Min(#rows,last_row);
print "TARGET_2248_TOR2244_ROWS101_500";
print "CONFIG","bank_file",bank_file,"bank_rows",#rows,
      "first",first,"last",last,"reduction_prime_count",reduction_prime_count;

base_curves:=[]; base_ai:=[];
candidate_curves:=[]; candidate_ai:=[];
source_curves:=[]; source_ai:=[]; source_fs:=[]; source_hs:=[];
source_rows:=[]; source_edges:=[]; source_cert_p:=[];
source_cert_chi:=[]; source_cert_degrees:=[];
source_status:=[]; source_cert_witness:=[]; source_product_edge:=[];

rows_seen:=0; rows_completed:=0; base_duplicates:=0; base_errors:=0;
reverse_edges:=0; reverse_jacobians:=0; reverse_one_split:=0;
reverse_gate_pass_raw:=0; reverse_candidate_duplicates:=0;
reverse_exact_tests:=0; reverse_exact_248:=0; reverse_cert_tests:=0;
reverse_simple_248:=0; reverse_product_248:=0; reverse_unresolved_248:=0;
source_duplicates:=0;

for ri in [first..last] do
    rows_seen+:=1; t:=rows[ri]; a,b,c,d:=Explode(t);
    try
        f0:=x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
        C0,f0n,h0n:=NormalizeCurve(HyperellipticCurve(f0));
        ai0:=AbsoluteInvariants(C0);
        oldbase:=ExistingCurve(C0,ai0,base_curves,base_ai);
        if oldbase ne 0 then
            base_duplicates+:=1;
            print "BASE_DUPLICATE","row",ri,"tuple",t,"existing",oldbase;
            rows_completed+:=1;
            continue;
        end if;
        Append(~base_curves,C0); Append(~base_ai,ai0);
        Rs:=RichelotIsogenousSurfaces(Jacobian(C0));
        print "BASE","row",ri,"tuple",t,"base_index",#base_curves,
              "richelot_count",#Rs;
        for j in [1..#Rs] do
            reverse_edges+:=1;
            if Type(Rs[j]) ne JacHyp then continue; end if;
            reverse_jacobians+:=1;
            D,fD,hD:=NormalizeCurve(Curve(Rs[j]));
            if not OneSplitPattern(D) then continue; end if;
            reverse_one_split+:=1;
            gate,g,redrows:=TwoAdicReductionGate(D,6);
            if not gate then continue; end if;
            reverse_gate_pass_raw+:=1;
            aiD:=AbsoluteInvariants(D);
            old:=ExistingCurve(D,aiD,candidate_curves,candidate_ai);
            if old ne 0 then
                reverse_candidate_duplicates+:=1;
                print "SOURCE_CANDIDATE_DUPLICATE","row",ri,"edge",j,
                      "existing_candidate",old,"reduction_gcd",g;
                continue;
            end if;
            Append(~candidate_curves,D); Append(~candidate_ai,aiD);
            reverse_exact_tests+:=1;
            G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
            print "SOURCE_CANDIDATE","row",ri,"tuple",t,"edge",j,
                  "candidate_index",#candidate_curves,
                  "factor_degrees",FactorDegrees(D),"reduction_gcd",g,
                  "reduction_rows",redrows,"torsion",inv,"order",#G,
                  "f",fD,"h",hD;
            if inv ne [2,4,8] then continue; end if;
            reverse_exact_248+:=1; reverse_cert_tests+:=1;
            simple,pcert,chi,degrees:=StrictRootPowerCertificate(D);
            witness:="source"; product:=false; product_edge:=0;
            if not simple then
                simple0,p0,chi0,degrees0:=StrictRootPowerCertificate(C0);
                if simple0 then
                    simple:=true; pcert:=p0; chi:=chi0; degrees:=degrees0;
                    witness:="reverse_base";
                else
                    product,product_edge:=RationalProductCertificate(D);
                    witness:=product select "rational_SetCart_quotient" else "none";
                end if;
            end if;
            status:=simple select "certified_simple" else
                    (product select "proven_product" else "unresolved");
            print "SOURCE_248_CERT","row",ri,"edge",j,
                  "status",status,"witness",witness,
                  "simple",BoolString(simple),"p",pcert,"chi",chi,
                  "degrees",degrees,"product",BoolString(product),
                  "product_edge",product_edge,"absolute_G2",aiD;
            if simple then reverse_simple_248+:=1;
            elif product then reverse_product_248+:=1;
            else reverse_unresolved_248+:=1;
            end if;
            oldsource:=ExistingCurve(D,aiD,source_curves,source_ai);
            if oldsource ne 0 then
                source_duplicates+:=1;
                print "SOURCE_248_DUPLICATE","row",ri,"edge",j,
                      "existing_source",oldsource;
                continue;
            end if;
            Append(~source_curves,D); Append(~source_ai,aiD);
            Append(~source_fs,fD); Append(~source_hs,hD);
            Append(~source_rows,ri); Append(~source_edges,j);
            Append(~source_cert_p,pcert); Append(~source_cert_chi,chi);
            Append(~source_cert_degrees,degrees);
            Append(~source_status,status); Append(~source_cert_witness,witness);
            Append(~source_product_edge,product_edge);
            print "SOURCE_ACCEPT","source",#source_curves,"row",ri,
                  "tuple",t,"edge",j,"f",fD,"h",hD,
                  "absolute_G2",aiD,"status",status,"witness",witness,
                  "cert_p",pcert,"cert_chi",chi,
                  "cert_degrees",degrees,"product_edge",product_edge;
        end for;
        rows_completed+:=1;
    catch e
        base_errors+:=1;
        print "BASE_ERROR","row",ri,"tuple",t,e`Object;
    end try;
    if rows_seen mod 10 eq 0 then
        print "PROGRESS","rows_seen",rows_seen,"rows_completed",rows_completed,
              "current_row",ri,"edges",reverse_edges,"one_split",reverse_one_split,
              "gate_pass",reverse_gate_pass_raw,"exact_tests",reverse_exact_tests,
              "exact_248",reverse_exact_248,"simple_sources",#source_curves;
    end if;
end for;

print "REVERSE_SUMMARY","requested_first",first,"requested_last",last,
      "rows_seen",rows_seen,"rows_completed",rows_completed,
      "base_unique",#base_curves,"base_duplicates",base_duplicates,
      "base_errors",base_errors,"edges",reverse_edges,
      "jacobians",reverse_jacobians,"one_split",reverse_one_split,
      "gate_pass_raw",reverse_gate_pass_raw,
      "candidate_duplicates",reverse_candidate_duplicates,
      "unique_exact_tests",reverse_exact_tests,"exact_248",reverse_exact_248,
      "cert_tests",reverse_cert_tests,"certified_simple_248",reverse_simple_248,
      "proven_product_248",reverse_product_248,
      "unresolved_248",reverse_unresolved_248,
      "source_duplicates",source_duplicates,
      "accepted_unique_exact_sources",#source_curves;

forward_edges:=0; forward_jacobians:=0; forward_full:=0;
forward_gate_pass:=0; forward_exact_tests:=0; forward_errors:=0;
target_hits:=0; neighbor_curves:=[]; neighbor_ai:=[];
neighbor_duplicates:=0;

for si in [1..#source_curves] do
    C:=source_curves[si];
    print "SOURCE_FORWARD","source",si,"row",source_rows[si],
          "origin_edge",source_edges[si],"f",source_fs[si],"h",source_hs[si],
          "status",source_status[si],"witness",source_cert_witness[si],
          "product_edge",source_product_edge[si],
          "cert_p",source_cert_p[si],"cert_chi",source_cert_chi[si],
          "cert_degrees",source_cert_degrees[si];
    try
        Rs:=RichelotIsogenousSurfaces(Jacobian(C));
        print "SOURCE_FORWARD_RICHELOT_COUNT","source",si,"count",#Rs;
        for j in [1..#Rs] do
            forward_edges+:=1;
            if Type(Rs[j]) ne JacHyp then
                print "FORWARD_EDGE","source",si,"edge",j,
                      "type",Type(Rs[j]);
                continue;
            end if;
            forward_jacobians+:=1;
            D,fD,hD:=NormalizeCurve(Curve(Rs[j])); ds:=FactorDegrees(D);
            full:=FullWeierstrassPattern(D); aiD:=AbsoluteInvariants(D);
            old:=ExistingCurve(D,aiD,neighbor_curves,neighbor_ai);
            if old ne 0 then neighbor_duplicates+:=1;
            else Append(~neighbor_curves,D); Append(~neighbor_ai,aiD); end if;
            print "FORWARD_EDGE","source",si,"edge",j,"type JacHyp",
                  "factor_degrees",ds,"full_weierstrass",BoolString(full),
                  "duplicate_neighbor",old;
            if not full then continue; end if;
            forward_full+:=1;
            gate,g,redrows:=TwoAdicReductionGate(D,7);
            if gate then forward_gate_pass+:=1; end if;
            // Exact testing is unconditional for every full-Weierstrass edge.
            forward_exact_tests+:=1;
            G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
            print "FORWARD_FULL_EXACT","source",si,"edge",j,
                  "gate_pass",BoolString(gate),"reduction_gcd",g,
                  "reduction_rows",redrows,"torsion",inv,"order",#G,
                  "f",fD,"h",hD,"absolute_G2",aiD;
            if inv eq [2,2,4,8] then
                target_hits+:=1;
                simpleT,pT,chiT,dsT:=StrictRootPowerCertificate(D);
                target_status:=source_status[si];
                target_witness:=Sprintf("isogenous_to_source_%o",si);
                target_product_edge:=0;
                if target_status eq "unresolved" then
                    if simpleT then
                        target_status:="certified_simple";
                        target_witness:="target_root_power";
                    else
                        prodT,target_product_edge:=RationalProductCertificate(D);
                        if prodT then
                            target_status:="proven_product";
                            target_witness:="target_rational_SetCart_quotient";
                        end if;
                    end if;
                end if;
                print "TARGET_2248_HIT","hit",target_hits,"source",si,
                      "source_row",source_rows[si],"source_edge",source_edges[si],
                      "forward_edge",j,"target_status",target_status,
                      "target_witness",target_witness,
                      "target_product_edge",target_product_edge,
                      "target_root_power",BoolString(simpleT),
                      "target_p",pT,"target_chi",chiT,"target_degrees",dsT,
                      "target_f",fD,"target_h",hD;
                PrintKernelIdentification(C,D,j);
            end if;
        end for;
    catch e
        forward_errors+:=1;
        print "SOURCE_FORWARD_ERROR","source",si,e`Object;
    end try;
end for;

print "FORWARD_SUMMARY","sources",#source_curves,"edges",forward_edges,
      "jacobians",forward_jacobians,"full_weierstrass",forward_full,
      "gate_pass",forward_gate_pass,"exact_tests",forward_exact_tests,
      "errors",forward_errors,"neighbor_duplicates",neighbor_duplicates,
      "unique_neighbors",#neighbor_curves,"target_hits",target_hits;
print "TARGET_2248_TOR2244_ROWS101_500_DONE";
quit;
