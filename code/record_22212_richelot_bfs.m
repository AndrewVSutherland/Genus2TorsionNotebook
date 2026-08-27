//////////////////////////////////////////////////////////////////////
// Rational Richelot graph of the order-96 [2,2,2,12] record.
//
// The collaborator's model is
//
//   y^2 + (1+x^2)y = f0(x).
//
// Completing the square gives Y^2 = (1+x^2)^2+4*f0.  This script
// factors that sextic, enumerates its 15 partitions into three pairs,
// constructs and identifies all 15 Richelot codomains, computes exact
// rational torsion and absolute G2 invariants, and then performs a
// rational Richelot BFS (default depth 2), deduplicating over Q.
//
// It also audits rational halves of the 15 nonzero source J[2]
// classes and the abstract exact torsion group.
//
// Typical run from the repository root:
//
//   magma -b code/record_22212_richelot_bfs.m
//
// Options:
//
//   max_depth:=1       only the 15 immediate codomains
//   write_log:=false   do not write the default results log
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned max_depth then
    max_depth := 2;
elif Type(max_depth) eq MonStgElt then
    max_depth := StringToInteger(max_depth);
end if;

if not assigned write_log then
    write_log := true;
elif Type(write_log) eq MonStgElt then
    write_log := write_log in {"true", "True", "1", "yes"};
end if;

if not assigned log_file then
    log_file := "results/record_22212_richelot_bfs.log";
end if;

if write_log then
    SetLogFile(log_file : Overwrite := true);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function InvOrder(inv)
    return #inv eq 0 select 1 else &*inv;
end function;

function InvExponent(inv)
    return #inv eq 0 select 1 else inv[#inv];
end function;

function IntegralSquareScale(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f), d;
end function;

function NormalizeCurve(C)
    f,h := HyperellipticPolynomials(C);
    // Complete the square if a generalized model occurs.
    F := h eq 0 select P!f else P!(h^2+4*f);
    FI,d := IntegralSquareScale(F);
    D := SimplifiedModel(HyperellipticCurve(FI));
    fD,hD := HyperellipticPolynomials(D);
    return D,P!fD,P!hD;
end function;

function TorsionData(C)
    D,fD,hD := NormalizeCurve(C);
    A,mp := TorsionSubgroup(Jacobian(D));
    return D,fD,hD,Invariants(A),A,mp;
end function;

function AbsoluteInvariants(C)
    return [Q!a : a in G2Invariants(C)];
end function;

function Bracket(A,B)
    return Derivative(A)*B-A*Derivative(B);
end function;

function CoefficientMatrixDeterminant(A,B,C)
    M := Matrix(Q,3,3,
        [Coefficient(A,i) : i in [0..2]] cat
        [Coefficient(B,i) : i in [0..2]] cat
        [Coefficient(C,i) : i in [0..2]]);
    return Determinant(M);
end function;

function AllPairings6()
    // Fix 1 in the first pair.  Pair it with j, then fix the least
    // remaining index in the second pair.  This gives 5*3=15 exactly.
    out := [];
    for j in [2..6] do
        rem := [k : k in [2..6] | k ne j];
        a := rem[1];
        for pos in [2..4] do
            b := rem[pos];
            rest := [k : k in rem | k ne a and k ne b];
            Append(~out,[<1,j>,<a,b>,<rest[1],rest[2]>]);
        end for;
    end for;
    assert #out eq 15;
    return out;
end function;

function PairQuadratic(roots,ij)
    return (x-roots[ij[1]])*(x-roots[ij[2]]);
end function;

function SameQCurve(C,D)
    try
        ok,mp := IsIsomorphic(C,D);
        return ok;
    catch e
        return false;
    end try;
end function;

function ExistingNode(C,absinv,node_curves,node_absinv)
    // Absolute invariants are the cheap geometric prefilter.  Retain
    // distinct Q-twists unless Magma proves them Q-isomorphic.
    for i in [1..#node_curves] do
        if node_absinv[i] ne absinv then
            continue;
        end if;
        if SameQCurve(C,node_curves[i]) then
            return i;
        end if;
    end for;
    return 0;
end function;

function TargetName(inv)
    targets := [
        <[4,12],"4_12">,
        <[2,4,12],"2_4_12">,
        <[2,24],"2_24">,
        <[2,2,24],"2_2_24">,
        <[2,2,2,24],"2_2_2_24">,
        <[2,2,4,12],"2_2_4_12">
    ];
    for rec in targets do
        if inv eq rec[1] then
            return rec[2];
        end if;
    end for;
    return "";
end function;

procedure PrintNode(tag,label,depth,parent,edge,C,f,h,inv,absinv)
    print tag,label,"depth",depth,"parent",parent,"edge",edge,
          "torsion",inv,"order",InvOrder(inv),
          "exponent",InvExponent(inv),
          "factor_degrees",Sort([Degree(t[1]) : t in Factorization(f)]);
    print tag cat "_CURVE",label,"f",f,"h",h;
    print tag cat "_ABSOLUTE_G2",label,absinv;
    target := TargetName(inv);
    if target ne "" then
        print "TARGET_HIT",target,label,"depth",depth,"torsion",inv,
              "curve_f",f,"curve_h",h;
    end if;
end procedure;

//////////////////////////////////////////////////////////////////////
// Source and completed-square factorization.
//////////////////////////////////////////////////////////////////////

f0 :=
    830742747091037849
    + 32014154874551031*x
    - 530648977741620*x^2
    - 15854483576121*x^3
    + 150572203590*x^4
    + 737595570*x^5
    + 756900*x^6;
h0 := 1+x^2;
F := h0^2+4*f0;

f_reported :=
    3027600*x^6 + 2950382280*x^5 + 602288814361*x^4
    - 63417934304484*x^3 - 2122595910966478*x^2
    + 128056619498204124*x + 3322970988364151397;

print "RECORD_22212_RICHELOT_BFS";
print "CONFIG","max_depth",max_depth,"write_log",write_log,
      "log_file",log_file;
print "GENERALIZED_F",f0;
print "GENERALIZED_H",h0;
print "COMPLETED_SQUARE_H2_PLUS_4F",F;
print "COMPLETED_SQUARE_MATCH_REPORTED",F eq f_reported;
assert F eq f_reported;
assert Discriminant(F) ne 0;

fac := Factorization(F);
print "COMPLETED_SQUARE_FACTORIZATION",fac;
assert #fac eq 6 and &and[t[2] eq 1 and Degree(t[1]) eq 1 : t in fac];
roots := Sort([-Coefficient(t[1],0)/Coefficient(t[1],1) : t in fac]);
print "SORTED_WEIERSTRASS_ROOTS",roots;
assert &*[x-r : r in roots] * LeadingCoefficient(F) eq F;

Csource,fs,hs,source_inv,Tsource,torsion_map :=
    TorsionData(HyperellipticCurve(F));
source_absinv := AbsoluteInvariants(Csource);
PrintNode("SOURCE","source",0,0,"input",Csource,fs,hs,
          source_inv,source_absinv);
assert source_inv eq [2,2,2,12];

//////////////////////////////////////////////////////////////////////
// Divisibility by 2 on the source.
//////////////////////////////////////////////////////////////////////

print "SOURCE_TORSION_GENERATOR_HALVING_START";
for i in [1..Ngens(Tsource)] do
    a := Tsource.i;
    Pa := torsion_map(a);
    divisible,half := IsDivisibleBy(Pa,2);
    print "SOURCE_TORSION_GENERATOR_HALVING",i,
          "abstract_order",Order(a),"divisible",divisible;
    if divisible then
        print "SOURCE_TORSION_GENERATOR_HALF",i,"half_order",Order(half),
              "check",2*half eq Pa;
    end if;
end for;

twice_set := Seqset([2*a : a in Tsource]);
print "SOURCE_TORSION_ABSTRACT_DOUBLE_IMAGE_SIZE",#twice_set;
for a in Tsource do
    if a in twice_set then
        halves := [b : b in Tsource | 2*b eq a];
        print "SOURCE_TORSION_ABSTRACT_DIVISIBLE",
              "coordinates",Eltseq(a),"order",Order(a),
              "number_of_halves",#halves,
              "one_half_coordinates",Eltseq(halves[1]),
              "one_half_order",Order(halves[1]);
    end if;
end for;

Jsource := Jacobian(Csource);
pairings := AllPairings6();
two_divisible_count := 0;
two_class_pairs := &cat[[<i,j> : j in [i+1..6]] : i in [1..5]];
assert #two_class_pairs eq 15;
print "SOURCE_J2_HALVING_START","class_count",#two_class_pairs;
for k in [1..#two_class_pairs] do
    // The 15 nonzero J[2] classes are indexed by unordered pairs of
    // Weierstrass roots (with complementary pairs giving the same
    // divisor class only when four roots, rather than six, are used).
    ij := two_class_pairs[k];
    u := PairQuadratic(roots,ij);
    P2 := Jsource![u,Q!0];
    assert P2 ne Jsource!0 and 2*P2 eq Jsource!0;
    divisible,half := IsDivisibleBy(P2,2);
    print "SOURCE_J2_HALVING_CLASS",k,"pair",ij,
          "roots",[roots[ij[1]],roots[ij[2]]],
          "u",u,"divisible",divisible;
    if divisible then
        two_divisible_count +:= 1;
        print "SOURCE_J2_HALVING_HIT",k,"pair",ij,
              "half_order",Order(half),"check",2*half eq P2;
    end if;
end for;
print "SOURCE_J2_HALVING_END","divisible_classes",two_divisible_count;
assert two_divisible_count eq 1;

//////////////////////////////////////////////////////////////////////
// Enumerate and identify all 15 source pairings.
//////////////////////////////////////////////////////////////////////

source_richelots := RichelotIsogenousSurfaces(Jsource);
print "SOURCE_RICHELOT_BUILTIN_COUNT",#source_richelots;
assert #source_richelots eq 15;
assert &and[Type(A) eq JacHyp : A in source_richelots];

builtin_curves := [];
builtin_f := [];
builtin_h := [];
builtin_absinv := [];
for i in [1..#source_richelots] do
    D,fD,hD := NormalizeCurve(Curve(source_richelots[i]));
    Append(~builtin_curves,D);
    Append(~builtin_f,fD);
    Append(~builtin_h,hD);
    Append(~builtin_absinv,AbsoluteInvariants(D));
end for;

node_curves := [Csource];
node_f := [fs];
node_h := [hs];
node_absinv := [source_absinv];
node_torsion := [source_inv];
node_depth := [0];
node_parent := [0];
node_edge := ["input"];
node_label := ["source"];

immediate_node_indices := [];
used_builtin := {};

for k in [1..#pairings] do
    pairing := pairings[k];
    A := PairQuadratic(roots,pairing[1]);
    B := PairQuadratic(roots,pairing[2]);
    C := PairQuadratic(roots,pairing[3]);
    Delta := CoefficientMatrixDeterminant(A,B,C);
    g0 := Delta*Bracket(B,C)*Bracket(C,A)*Bracket(A,B);
    smooth := Delta ne 0 and Degree(g0) in {5,6} and Discriminant(g0) ne 0;
    print "SOURCE_PAIRING",k,"pairs",pairing,
          "root_pairs",[
              [roots[ij[1]],roots[ij[2]]] : ij in pairing
          ],"quadratics",[A,B,C],"Delta",Delta,"smooth",smooth;
    assert smooth;

    Cplus,fp,hp := NormalizeCurve(HyperellipticCurve(g0));
    Cminus,fm,hm := NormalizeCurve(HyperellipticCurve(-g0));
    pair_absinv := AbsoluteInvariants(Cplus);
    assert pair_absinv eq AbsoluteInvariants(Cminus);

    matches := [i : i in [1..#builtin_curves] |
        builtin_absinv[i] eq pair_absinv];
    assert #matches eq 1;
    bi := matches[1];
    assert bi notin used_builtin;
    Include(~used_builtin,bi);

    plus_iso := SameQCurve(Cplus,builtin_curves[bi]);
    minus_iso := SameQCurve(Cminus,builtin_curves[bi]);
    assert plus_iso xor minus_iso;
    if plus_iso then
        D := Cplus; fD := fp; hD := hp; formula_sign := 1;
    else
        D := Cminus; fD := fm; hD := hm; formula_sign := -1;
    end if;

    DT,fT,hT,invT,AT,mpT := TorsionData(D);
    absT := AbsoluteInvariants(DT);
    label := Sprintf("pairing_%o",k);
    edge := Sprintf("source_pairing_%o_builtin_%o_sign_%o",k,bi,formula_sign);
    print "SOURCE_PAIRING_IDENTIFICATION",k,"builtin_index",bi,
          "formula_sign",formula_sign,"manual_g0",g0,
          "plus_iso",plus_iso,"minus_iso",minus_iso;

    existing := ExistingNode(DT,absT,node_curves,node_absinv);
    if existing eq 0 then
        Append(~node_curves,DT);
        Append(~node_f,fT);
        Append(~node_h,hT);
        Append(~node_absinv,absT);
        Append(~node_torsion,invT);
        Append(~node_depth,1);
        Append(~node_parent,1);
        Append(~node_edge,edge);
        Append(~node_label,label);
        ni := #node_curves;
        Append(~immediate_node_indices,ni);
        PrintNode("IMMEDIATE",label,1,1,edge,DT,fT,hT,invT,absT);
    else
        Append(~immediate_node_indices,existing);
        print "IMMEDIATE_DUPLICATE",label,"existing_node",existing,
              "torsion",invT,"absolute_G2",absT;
        assert node_torsion[existing] eq invT;
    end if;
end for;
assert #used_builtin eq 15;
print "SOURCE_PAIRINGS_ALL_IDENTIFIED",#used_builtin;
print "IMMEDIATE_UNIQUE_Q_CURVES",#Seqset(immediate_node_indices);

//////////////////////////////////////////////////////////////////////
// Breadth-first traversal beyond the immediate layer.
//////////////////////////////////////////////////////////////////////

edges_total := 15;
duplicate_edges := 15-#Seqset(immediate_node_indices);
nonjac_edges := 0;

if max_depth ge 2 then
    cursor := 2;
    while cursor le #node_curves do
        d := node_depth[cursor];
        if d ge max_depth then
            cursor +:= 1;
            continue;
        end if;

        print "BFS_EXPAND",cursor,node_label[cursor],"depth",d,
              "torsion",node_torsion[cursor];
        Rs := RichelotIsogenousSurfaces(Jacobian(node_curves[cursor]));
        print "BFS_RICHELOT_COUNT",cursor,#Rs;
        for j in [1..#Rs] do
            edges_total +:= 1;
            edge := Sprintf("node_%o_richelot_%o",cursor,j);
            if Type(Rs[j]) ne JacHyp then
                nonjac_edges +:= 1;
                print "BFS_NONJAC_EDGE",cursor,j,"type",Type(Rs[j]),Rs[j];
                continue;
            end if;

            D,fD,hD := NormalizeCurve(Curve(Rs[j]));
            absD := AbsoluteInvariants(D);
            existing := ExistingNode(D,absD,node_curves,node_absinv);
            if existing ne 0 then
                duplicate_edges +:= 1;
                print "BFS_DUPLICATE_EDGE","parent",cursor,"edge_index",j,
                      "existing_node",existing,"absolute_G2",absD;
                continue;
            end if;

            DT,fT,hT,invT,AT,mpT := TorsionData(D);
            absT := AbsoluteInvariants(DT);
            // Simplification is a Q-isomorphism, so this must retain the key.
            assert absT eq absD;
            label := Sprintf("node_%o",#node_curves+1);
            Append(~node_curves,DT);
            Append(~node_f,fT);
            Append(~node_h,hT);
            Append(~node_absinv,absT);
            Append(~node_torsion,invT);
            Append(~node_depth,d+1);
            Append(~node_parent,cursor);
            Append(~node_edge,edge);
            Append(~node_label,label);
            PrintNode("BFS_NEW",label,d+1,cursor,edge,DT,fT,hT,invT,absT);
        end for;
        cursor +:= 1;
    end while;
end if;

//////////////////////////////////////////////////////////////////////
// Summaries.
//////////////////////////////////////////////////////////////////////

print "TORSION_TYPE_SUMMARY_START";
torsion_types := [];
for inv in node_torsion do
    if inv notin torsion_types then
        Append(~torsion_types,inv);
    end if;
end for;
for inv in torsion_types do
    ids := [i : i in [1..#node_torsion] | node_torsion[i] eq inv];
    print "TORSION_TYPE",inv,"order",InvOrder(inv),
          "count",#ids,"nodes",ids;
end for;
print "TORSION_TYPE_SUMMARY_END";

print "IMMEDIATE_TORSION_TYPE_SUMMARY_START";
immediate_types := [];
for i in Seqset(immediate_node_indices) do
    inv := node_torsion[i];
    if inv notin immediate_types then Append(~immediate_types,inv); end if;
end for;
for inv in immediate_types do
    ids := [i : i in Seqset(immediate_node_indices) |
        node_torsion[i] eq inv];
    print "IMMEDIATE_TORSION_TYPE",inv,"order",InvOrder(inv),
          "count",#ids,"nodes",ids;
end for;
print "IMMEDIATE_TORSION_TYPE_SUMMARY_END";

target_hits := 0;
for i in [1..#node_curves] do
    target := TargetName(node_torsion[i]);
    if target ne "" then
        target_hits +:= 1;
        print "FINAL_TARGET_HIT",target,"node",i,node_label[i],
              "depth",node_depth[i],"torsion",node_torsion[i],
              "curve_f",node_f[i],"curve_h",node_h[i],
              "absolute_G2",node_absinv[i];
    end if;
end for;

print "BFS_SUMMARY","max_depth",max_depth,"nodes",#node_curves,
      "edges",edges_total,"duplicate_edges",duplicate_edges,
      "nonjac_edges",nonjac_edges,"target_hits",target_hits;
print "GEOMETRIC_SIMPLICITY_INHERITANCE",
      "Every Jacobian node is Q-isogenous to the source by a chain of Richelot isogenies; geometric simplicity is isogeny invariant.";
print "SOURCE_SIMPLICITY_CERTIFICATE",
      "code/verify_record_22212_order96.m",
      "with independent root-power Frobenius witnesses at p=37 and p=73.";
print "DONE_RECORD_22212_RICHELOT_BFS";

if write_log then
    UnsetLogFile();
end if;
quit;
