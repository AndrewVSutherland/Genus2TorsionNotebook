//////////////////////////////////////////////////////////////////////
// Rational Richelot neighborhood of the certified simple [2,4,8] seed.
// Target: a Q-isogenous Jacobian with exact torsion [2,2,4,8].
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned max_depth then
    max_depth := 2;
elif Type(max_depth) eq MonStgElt then
    max_depth := StringToInteger(max_depth);
end if;

if not assigned max_nodes then
    max_nodes := 100;
elif Type(max_nodes) eq MonStgElt then
    max_nodes := StringToInteger(max_nodes);
end if;

if not assigned write_log then
    write_log := true;
elif Type(write_log) eq MonStgElt then
    write_log := write_log in {"true","True","1","yes"};
end if;

if not assigned log_file then
    log_file := "results/target_2248_richelot_seed.log";
end if;
if write_log then SetLogFile(log_file : Overwrite := true); end if;

Q := Rationals(); Z := Integers(); P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f);
end function;

function NormalizeCurve(C)
    f,h := HyperellipticPolynomials(C);
    F := h eq 0 select P!f else P!(h^2+4*f);
    FI := IntegralSquareScale(F);
    D := SimplifiedModel(HyperellipticCurve(FI));
    fD,hD := HyperellipticPolynomials(D);
    return D,P!fD,P!hD;
end function;

function TorsionData(C)
    D,f,h := NormalizeCurve(C);
    G,mp := TorsionSubgroup(Jacobian(D));
    return D,f,h,Invariants(G);
end function;

function AbsoluteInvariants(C)
    return [Q!a : a in G2Invariants(C)];
end function;

function SameQCurve(C,D)
    try
        ok,mp := IsIsomorphic(C,D);
        return ok;
    catch er
        return false;
    end try;
end function;

function ExistingNode(C,ai,curves,ais)
    for i in [1..#curves] do
        if ais[i] eq ai and SameQCurve(C,curves[i]) then return i; end if;
    end for;
    return 0;
end function;

function RootPowerCertificate(f)
    C := HyperellipticCurve(f);
    for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
               67,71,73,79,83,89,97,101,103,107,109,113,127,131] do
        try
            fp := ChangeRing(f,GF(pp));
            if Degree(fp) ne Degree(f) or Discriminant(fp) eq 0 then continue; end if;
            Lp := LPolynomial(ChangeRing(C,GF(pp)));
            cs := [Z!Coefficient(Lp,i) : i in [0..4]];
            S<T> := PolynomialRing(Q);
            chi := &+[Q!cs[5-i]*T^i : i in [0..4]];
            if not IsIrreducible(chi) then continue; end if;
            K<pi> := NumberField(chi);
            ds := [Degree(MinimalPolynomial(pi^n)) : n in [2..12]];
            if &and[d eq 4 : d in ds] then return true,pp,chi,ds; end if;
        catch er
            continue;
        end try;
    end for;
    return false,0,P!0,[];
end function;

procedure PrintNode(i,depth,parent,edge,f,h,inv)
    print "NODE",i,"depth",depth,"parent",parent,"edge",edge,
          "torsion",inv,"order",&*inv,
          "factor_degrees",Sort([Degree(a[1]) : a in Factorization(f)]);
    print "NODE_CURVE",i,"f",f,"h",h;
end procedure;

// Certified simple exact [2,4,8] seed from notes/m18_m14_halving.md.
fseed :=
    7061463847622250*x^5
    +104632219276049025*x^4
    +135735215960638800*x^3
    +188573481843278400*x^2
    +51200550567936000*x;

C0,f0,h0,inv0 := TorsionData(HyperellipticCurve(fseed));
assert inv0 eq [2,4,8];
ok0,p0,chi0,ds0 := RootPowerCertificate(f0);

print "TARGET_2248_RICHELOT_SEED";
print "CONFIG","max_depth",max_depth,"max_nodes",max_nodes;
print "SOURCE_ROOT_POWER","p",p0,"chi",chi0,"degrees",ds0;
if not ok0 then
    print "SOURCE_ROOT_POWER_STRICT_NOT_FOUND_IN_PRIME_LIST";
    for pp in [47] do
        fp := ChangeRing(f0,GF(pp));
        print "SOURCE_DIAGNOSTIC","p",pp,"degree",Degree(fp),
              "discriminant_zero",Discriminant(fp) eq 0;
        if Degree(fp) eq Degree(f0) and Discriminant(fp) ne 0 then
            Lp := LPolynomial(ChangeRing(HyperellipticCurve(f0),GF(pp)));
            cs := [Z!Coefficient(Lp,i) : i in [0..4]];
            Sd<Td> := PolynomialRing(Q);
            chid := &+[Q!cs[5-i]*Td^i : i in [0..4]];
            Kd<pid> := NumberField(chid);
            print "SOURCE_DIAGNOSTIC_LP",Lp,"chi",chid,
                  "degrees",[Degree(MinimalPolynomial(pid^n)) : n in [2..12]];
        end if;
    end for;
end if;

curves := [C0]; fs := [f0]; hs := [h0]; invs := [inv0];
ais := [AbsoluteInvariants(C0)];
depths := [0]; parents := [0]; edges := ["source"];
PrintNode(1,0,0,"source",f0,h0,inv0);

cursor := 1;
edge_count := 0;
target_count := 0;
product_count := 0;

while cursor le #curves and #curves lt max_nodes do
    if depths[cursor] ge max_depth then cursor +:= 1; continue; end if;
    try
        Rs := RichelotIsogenousSurfaces(Jacobian(curves[cursor]));
        print "EXPAND",cursor,"depth",depths[cursor],"richelot_count",#Rs;
        for j in [1..#Rs] do
            edge_count +:= 1;
            edge := Sprintf("node_%o_richelot_%o",cursor,j);
            if Type(Rs[j]) ne JacHyp then
                product_count +:= 1;
                print "RICHELOT_NONJACOBIAN",edge,"type",Type(Rs[j]);
                continue;
            end if;
            D,fD,hD,invD := TorsionData(Curve(Rs[j]));
            aiD := AbsoluteInvariants(D);
            old := ExistingNode(D,aiD,curves,ais);
            if old ne 0 then
                print "RICHELOT_OLD",edge,"node",old,"torsion",invD;
                continue;
            end if;
            Append(~curves,D); Append(~fs,fD); Append(~hs,hD);
            Append(~invs,invD); Append(~ais,aiD);
            Append(~depths,depths[cursor]+1); Append(~parents,cursor);
            Append(~edges,edge);
            ni := #curves;
            PrintNode(ni,depths[ni],cursor,edge,fD,hD,invD);
            if invD eq [2,2,4,8] then
                target_count +:= 1;
                oks,ps,chis,dss := RootPowerCertificate(fD);
                print "TARGET_HIT",ni,"torsion",invD,
                      "root_power",oks,"p",ps,"chi",chis,"degrees",dss;
            end if;
            if #curves ge max_nodes then break; end if;
        end for;
    catch er
        print "EXPAND_ERROR",cursor,er`Object;
    end try;
    cursor +:= 1;
end while;

max_depth_reached,max_depth_index := Max(depths);
print "SUMMARY","nodes",#curves,"edges_examined",edge_count,
      "products",product_count,"targets",target_count,
      "max_depth_reached",max_depth_reached;
print "TORSION_DISTRIBUTION";
keys := Sort(Setseq(Seqset([Sprint(v) : v in invs])));
for key in keys do
    print key,#[i : i in [1..#invs] | Sprint(invs[i]) eq key];
end for;
print "ISOGENY_SIMPLICITY_NOTE all Jacobian nodes are geometrically simple because the source is and geometric simplicity is isogeny invariant";
print "TARGET_2248_RICHELOT_SEED_DONE";
quit;
