//////////////////////////////////////////////////////////////////////
// Focused audit of the two [2,4,8] source hits found at row 51,
// edges 13 and 14, of tor2244.txt.
//
// The script:
//   * reconstructs both sources from the recorded bank row;
//   * computes exact torsion and tests Q-isomorphism;
//   * gives a strict root-power simplicity certificate when available;
//   * identifies the reverse Richelot kernel from the fully split base;
//   * enumerates every rational forward Richelot neighbor;
//   * runs the full-J[2], v_2(#J(F_p)) and exact-torsion controls; and
//   * identifies every forward kernel by its quadratic splitting.
//
// Run from torsion_jac:
//   magma -b write_log:=true code/target_2248_source51_forward.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned bank_file then bank_file:="data/tor2244_bank.txt"; end if;
if not assigned write_log then write_log:=true;
elif Type(write_log) eq MonStgElt then write_log:=write_log in {"true","True","1","yes"}; end if;
if not assigned log_file then
    log_file:="results/target_2248_source51_forward.log";
end if;
if write_log then SetLogFile(log_file : Overwrite:=true); end if;

Q:=Rationals(); Z:=Integers(); P<x>:=PolynomialRing(Q);
prime_list:=[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
             67,71,73,79,83,89,97,101,103,107,109,113,127,131];

function BoolString(b) return b select "true" else "false"; end function;

function ReadTupleFile(filename)
    rows:=Split(Read(filename),"\n"); out:=[];
    for raw in rows do
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

function FactorDegrees(C)
    F:=CompletedPolynomial(C);
    return Sort(&cat[[Degree(z[1]):j in [1..z[2]]]:z in Factorization(F)]);
end function;

function FullWeierstrassPattern(C)
    F:=CompletedPolynomial(C); ds:=FactorDegrees(C);
    if Degree(F) eq 5 then return ds eq [1,1,1,1,1]; end if;
    return Degree(F) eq 6 and ds eq [1,1,1,1,1,1];
end function;

function AbsoluteInvariants(C)
    return [Q!a:a in G2Invariants(C)];
end function;

function SameQCurve(C,D)
    try ok,mp:=IsIsomorphic(C,D); return ok;
    catch e return false; end try;
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
            if used ge 2 then return true,g,rows; end if;
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

function MatchKernelTriple(Ctarget,A,B,D,leading_scale)
    delta:=CoefficientMatrixDeterminant(A,B,D);
    if delta eq 0 then return false,Q!0,delta,P!0; end if;
    g0:=delta*Bracket(B,D)*Bracket(D,A)*Bracket(A,B);
    if Discriminant(g0) eq 0 then return false,Q!0,delta,g0; end if;
    // Factorization() returns monic factors and suppresses lc(F).  The
    // built-in Richelot codomain can therefore differ from the bare bracket
    // formula by the source leading coefficient (and by sign).
    scales:=[Q!1,Q!-1,Q!leading_scale,Q!-leading_scale];
    for u in scales do
        Cu,fu,hu:=NormalizeCurve(HyperellipticCurve(u*g0));
        if SameQCurve(Cu,Ctarget) then return true,u,delta,g0; end if;
    end for;
    return false,Q!0,delta,g0;
end function;

// Identify one of the 15 kernels on a fully split genus-2 curve.
procedure PrintFullSplitKernelIdentification(Csource,Ctarget,builtin_index)
    F:=CompletedPolynomial(Csource); fac:=Factorization(F);
    lins:=[z[1]:z in fac|Degree(z[1]) eq 1];
    triples:=[]; labels:=[];
    if Degree(F) eq 5 and #lins eq 5 then
        for k in [1..5] do
            rest:=[i:i in [1..5]|i ne k];
            for pairing in Pairings4() do
                i1:=rest[pairing[1][1]]; i2:=rest[pairing[1][2]];
                i3:=rest[pairing[2][1]]; i4:=rest[pairing[2][2]];
                Append(~triples,<lins[i1]*lins[i2],lins[i3]*lins[i4],lins[k]>);
                Append(~labels,Sprintf("L%o*L%o | L%o*L%o | L%o*infinity",
                                      i1,i2,i3,i4,k));
            end for;
        end for;
    elif Degree(F) eq 6 and #lins eq 6 then
        matchings:=[
          [<1,2>,<3,4>,<5,6>],[<1,2>,<3,5>,<4,6>],[<1,2>,<3,6>,<4,5>],
          [<1,3>,<2,4>,<5,6>],[<1,3>,<2,5>,<4,6>],[<1,3>,<2,6>,<4,5>],
          [<1,4>,<2,3>,<5,6>],[<1,4>,<2,5>,<3,6>],[<1,4>,<2,6>,<3,5>],
          [<1,5>,<2,3>,<4,6>],[<1,5>,<2,4>,<3,6>],[<1,5>,<2,6>,<3,4>],
          [<1,6>,<2,3>,<4,5>],[<1,6>,<2,4>,<3,5>],[<1,6>,<2,5>,<3,4>]
        ];
        for m in matchings do
            Append(~triples,<lins[m[1][1]]*lins[m[1][2]],
                             lins[m[2][1]]*lins[m[2][2]],
                             lins[m[3][1]]*lins[m[3][2]]>);
            Append(~labels,Sprintf("L%o*L%o | L%o*L%o | L%o*L%o",
                m[1][1],m[1][2],m[2][1],m[2][2],m[3][1],m[3][2]));
        end for;
    else
        print "REVERSE_KERNEL_IDENTIFICATION_UNAVAILABLE","factorization",fac;
        return;
    end if;
    matches:=0;
    for k in [1..#triples] do
        A:=triples[k][1]; B:=triples[k][2]; D:=triples[k][3];
        ok,sgn,delta,g0:=MatchKernelTriple(Ctarget,A,B,D,LeadingCoefficient(F));
        if ok then
            matches+:=1;
            print "SOURCE_REVERSE_KERNEL_IDENTIFICATION","builtin_index",builtin_index,
                  "kernel_index",k,"label",labels[k],"factors",[A,B,D],
                  "delta",delta,"formula_scale",sgn,"manual_g0",g0;
        end if;
    end for;
    print "SOURCE_REVERSE_KERNEL_MATCH_COUNT",builtin_index,matches;
end procedure;

// Identify one of the three rational kernels on a one-split source.
procedure PrintOneSplitKernelIdentification(Csource,Ctarget,builtin_index)
    F:=CompletedPolynomial(Csource); fac:=Factorization(F);
    lins:=[z[1]:z in fac|Degree(z[1]) eq 1];
    quads:=[z[1]:z in fac|Degree(z[1]) eq 2];
    triples:=[]; labels:=[];
    if #quads ne 1 then
        print "FORWARD_KERNEL_IDENTIFICATION_UNAVAILABLE",builtin_index,
              "factorization",fac;
        return;
    end if;
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
        print "FORWARD_KERNEL_IDENTIFICATION_UNAVAILABLE",builtin_index,
              "factorization",fac;
        return;
    end if;
    matches:=0;
    for k in [1..#triples] do
        A:=triples[k][1]; B:=triples[k][2]; D:=triples[k][3];
        ok,sgn,delta,g0:=MatchKernelTriple(Ctarget,A,B,D,LeadingCoefficient(F));
        if ok then
            matches+:=1;
            print "FORWARD_KERNEL_IDENTIFICATION","builtin_index",builtin_index,
                  "kernel_index",k,"label",labels[k],"factors",[A,B,D],
                  "delta",delta,"formula_scale",sgn,"manual_g0",g0;
        end if;
    end for;
    print "FORWARD_KERNEL_MATCH_COUNT",builtin_index,matches;
end procedure;

print "TARGET_2248_SOURCE51_FOCUSED_FORWARD";
rows:=ReadTupleFile(bank_file); ri:=51; t:=rows[ri];
a,b,c,d:=Explode(t);
f0:=x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
C0,f0n,h0n:=NormalizeCurve(HyperellipticCurve(f0));
G0,mp0:=TorsionSubgroup(Jacobian(C0));
simple0,p0,chi0,degrees0:=StrictRootPowerCertificate(C0);
print "BASE","row",ri,"tuple",t,"f",f0n,"h",h0n,
      "torsion",Invariants(G0),"order",#G0,"factor_degrees",FactorDegrees(C0),
      "root_power",BoolString(simple0),"root_power_p",p0,
      "root_power_chi",chi0,"root_power_degrees",degrees0;

Rs0:=RichelotIsogenousSurfaces(Jacobian(C0));
print "BASE_RICHELOT_COUNT",#Rs0;
raw_sources:=[]; raw_fs:=[]; raw_hs:=[]; raw_edges:=[];
for j in [13,14] do
    if Type(Rs0[j]) ne JacHyp then
        print "SOURCE_EDGE_NOT_JACHYP",j,"type",Type(Rs0[j]);
        continue;
    end if;
    D,fD,hD:=NormalizeCurve(Curve(Rs0[j]));
    G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
    simple,pcert,chi,degrees:=StrictRootPowerCertificate(D);
    print "SOURCE_RAW","row",ri,"edge",j,"f",fD,"h",hD,
          "factor_degrees",FactorDegrees(D),"torsion",inv,"order",#G,
          "absolute_G2",AbsoluteInvariants(D),"root_power",BoolString(simple),
          "root_power_p",pcert,"root_power_chi",chi,"root_power_degrees",degrees;
    PrintFullSplitKernelIdentification(C0,D,j);
    Append(~raw_sources,D); Append(~raw_fs,fD); Append(~raw_hs,hD);
    Append(~raw_edges,j);
end for;

if #raw_sources eq 2 then
    print "SOURCE_PAIR_Q_ISOMORPHIC",BoolString(SameQCurve(raw_sources[1],raw_sources[2]));
end if;

sources:=[]; source_fs:=[]; source_hs:=[]; source_edges:=[];
for i in [1..#raw_sources] do
    old:=0;
    for k in [1..#sources] do
        if SameQCurve(raw_sources[i],sources[k]) then old:=k; break; end if;
    end for;
    if old ne 0 then
        print "SOURCE_Q_ISOMORPHISM_DUPLICATE","raw",i,"edge",raw_edges[i],
              "unique_source",old,"existing_edge",source_edges[old];
    else
        Append(~sources,raw_sources[i]); Append(~source_fs,raw_fs[i]);
        Append(~source_hs,raw_hs[i]); Append(~source_edges,raw_edges[i]);
        print "SOURCE_UNIQUE_ACCEPT",#sources,"raw",i,"edge",raw_edges[i];
    end if;
end for;
print "SOURCE_DEDUP_SUMMARY","raw",#raw_sources,"unique_Q_curves",#sources;

forward_edges:=0; forward_jacobians:=0; forward_full:=0;
forward_gate_pass:=0; forward_exact:=0; target_hits:=0;
for si in [1..#sources] do
    C:=sources[si];
    print "SOURCE_FORWARD",si,"reverse_edge",source_edges[si],
          "f",source_fs[si],"h",source_hs[si];
    Rs:=RichelotIsogenousSurfaces(Jacobian(C));
    print "SOURCE_FORWARD_RICHELOT_COUNT",si,#Rs;
    for j in [1..#Rs] do
        forward_edges+:=1;
        if Type(Rs[j]) ne JacHyp then
            print "FORWARD_EDGE",si,j,"type",Type(Rs[j]),"surface",Rs[j];
            if Type(Rs[j]) eq SetCart then
                print "SOURCE_NONSIMPLE_CERTIFICATE","source",si,"edge",j,
                      "reason","rational Richelot quotient is an elliptic-curve product",
                      "surface",Rs[j];
            end if;
            continue;
        end if;
        forward_jacobians+:=1;
        D,fD,hD:=NormalizeCurve(Curve(Rs[j]));
        ds:=FactorDegrees(D); full:=FullWeierstrassPattern(D);
        if full then forward_full+:=1; end if;
        gate,g,redrows:=TwoAdicReductionGate(D,7);
        if gate then forward_gate_pass+:=1; end if;
        // Exact torsion is deliberately run on every rational Jacobian
        // neighbor, stronger than the exact_all_full control.
        G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
        forward_exact+:=1;
        print "FORWARD_EXACT","source",si,"edge",j,"factor_degrees",ds,
              "full_weierstrass",BoolString(full),"v2_gate_7",BoolString(gate),
              "reduction_gcd",g,"reduction_rows",redrows,"torsion",inv,
              "order",#G,"f",fD,"h",hD,"absolute_G2",AbsoluteInvariants(D);
        PrintOneSplitKernelIdentification(C,D,j);
        if inv eq [2,2,4,8] then
            target_hits+:=1;
            simple,pT,chiT,degreesT:=StrictRootPowerCertificate(D);
            print "TARGET_2248_HIT","source",si,"edge",j,"target_f",fD,
                  "target_h",hD,"target_root_power",BoolString(simple),
                  "target_p",pT,"target_chi",chiT,
                  "target_root_power_degrees",degreesT;
        end if;
    end for;
end for;

print "FORWARD_SUMMARY","sources",#sources,"edges",forward_edges,
      "jacobians",forward_jacobians,"full_weierstrass",forward_full,
      "v2_gate_7_pass",forward_gate_pass,"exact_tests",forward_exact,
      "target_hits",target_hits;
print "TARGET_2248_SOURCE51_FOCUSED_FORWARD_DONE";
quit;
