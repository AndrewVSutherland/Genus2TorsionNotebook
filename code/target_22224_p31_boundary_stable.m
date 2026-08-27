//////////////////////////////////////////////////////////////////////
// Stable special-fibre test for the p=31 boundary of the q-square
// [2,2,2,24] cover.
//
// On a boundary point having exactly one double branch point, the
// normalization is an elliptic curve E and Pic^0 of the nodal fibre is an
// extension of E by the split or nonsplit one-dimensional torus.  This
// program tests divisibility by 2 in that COMPLETE generalized Jacobian,
// provided the marked Mumford divisor and the relevant tangent functions
// avoid the two points above the node.  Presentations meeting the node are
// reported separately: their answer depends on a first-order blow-up chart.
//
// Run:
//   magma -b code/target_22224_p31_boundary_stable.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned p then p := 31; end if;
if not assigned log_file then
    log_file := "results/target_22224_p31_boundary_stable.log";
end if;
SetLogFile(log_file : Overwrite:=true);

K := GF(p);  L := GF(p^2);  Z := Integers();
P<X> := PolynomialRing(K);

function MakeMonic(g)
    return g/LeadingCoefficient(g);
end function;

function DirectData(aa,bb,cc,dd,uu,tt,vv)
    f := X*(X+aa^2)*(X+bb^2)*(X+cc^2)*(X+dd^2);
    q := X^2+uu*X+tt^2;
    h := X^3+(1+3*uu)/2*X^2+vv*X+tt^3;
    return f,q,h;
end function;

// This is deliberately valid on a nodal fibre.  Only the CRT/contact chart
// is required to be open; smoothness of f is not required here.
function MarkedD12Data(aa,bb,cc,dd,uu,tt,vv)
    f,q,h := DirectData(aa,bb,cc,dd,uu,tt,vv);
    if h^2-f-q^3 ne 0 then return false,f,P!0,P!0,"contact"; end if;
    g4 := (X-aa*bb)*(X-cc*dd)-X*(aa+bb)*(cc+dd);
    L4 := (X-aa*bb)*(cc+dd)+(aa+bb)*(X-cc*dd);
    if Degree(g4) ne 2 then return false,f,P!0,P!0,"g4_degree"; end if;
    v3 := h mod q;  v4 := (-X*L4) mod g4;
    Rm := (q*X) mod g4;  Rn := q mod g4;  Rt := (v4-v3) mod g4;
    a0:=Coefficient(Rm,0); a1:=Coefficient(Rm,1);
    b0:=Coefficient(Rn,0); b1:=Coefficient(Rn,1);
    c0:=Coefficient(Rt,0); c1:=Coefficient(Rt,1);
    det := a0*b1-a1*b0;
    if det eq 0 then return false,f,P!0,P!0,"CRT_det"; end if;
    mc := (c0*b1-c1*b0)/det;
    nc := (a0*c1-a1*c0)/det;
    ell := v3+q*(mc*X+nc);
    if mc eq 0 or (ell^2-f) mod (q*g4) ne 0 then
        return false,f,P!0,P!0,"CRT_line";
    end if;
    Uraw := ExactQuotient(ell^2-f,q*g4);
    if Degree(Uraw) ne 2 then return false,f,P!0,P!0,"U_degree"; end if;
    U := MakeMonic(Uraw);  V := (-ell) mod U;
    if (V^2-f) mod U ne 0 then return false,f,P!0,P!0,"Mumford"; end if;
    return true,f,U,V,"open";
end function;

function RootsBoth(z)
    if z eq 0 then return [K!0]; end if;
    r := SquareRoot(z);
    return [r,-r];
end function;

function PointClass(J,E,pt)
    if pt[3] eq 0 then return J!0; end if;
    xx := pt[1]/pt[3];
    yy := pt[2]/pt[3]^2;
    return J![X-xx,yy];
end function;

// The CRT Mumford representative of D3+D4 can hit the node or cease to be
// an open chart even when D3 and D4 separately extend across the
// normalization.  This routine computes their ordinary Picard projection
// directly.  Failure of divisibility in Pic^0(E) is already a rigorous
// obstruction in the full generalized Jacobian.
function ProjectionD3D4(aa,bb,cc,dd,uu,tt,vv,node)
    f,q,h:=DirectData(aa,bb,cc,dd,uu,tt,vv);
    g4:=(X-aa*bb)*(X-cc*dd)-X*(aa+bb)*(cc+dd);
    L4:=(X-aa*bb)*(cc+dd)+(aa+bb)*(X-cc*dd);
    if Degree(g4) ne 2 then return false,"projection_g4",0,0; end if;
    if Degree(GCD(q,X-node)) gt 0 or Degree(GCD(g4,X-node)) gt 0 then
        return false,"projection_divisor_hits_node",0,0;
    end if;
    g:=ExactQuotient(f,(X-node)^2);
    if Degree(g) ne 3 or Discriminant(g) eq 0 then
        return false,"projection_not_one_node",0,0;
    end if;
    v3:=h mod q; v4:=(-X*L4) mod g4;
    V3:=(v3*InverseMod(X-node,q)) mod q;
    V4:=(v4*InverseMod(X-node,g4)) mod g4;
    U3:=MakeMonic(q); U4:=MakeMonic(g4);
    if (V3^2-g) mod U3 ne 0 or (V4^2-g) mod U4 ne 0 then
        return false,"projection_normalization_map",0,0;
    end if;
    E:=HyperellipticCurve(g); J:=Jacobian(E);
    D:=J![U3,V3]+J![U4,V4];
    nhalf:=0;
    for pt in Points(E) do
        if 2*PointClass(J,E,pt) eq D then nhalf+:=1; end if;
    end for;
    return true,"projection_complete",Order(D),nhalf;
end function;

// Returns whether the marked class is divisible by 2 in Pic^0 of the
// one-node special fibre.  This includes the torus squareclass, not merely
// divisibility of its elliptic projection.
function StableDivisible(f,U,V,node)
    if Degree(GCD(U,X-node)) gt 0 then
        return false,"marked_hits_node",0,0,0;
    end if;
    g := ExactQuotient(f,(X-node)^2);
    if Degree(g) ne 3 or Discriminant(g) eq 0 or Evaluate(g,node) eq 0 then
        return false,"not_one_node",0,0,0;
    end if;
    Vn := (V*InverseMod(X-node,U)) mod U;
    if (Vn^2-g) mod U ne 0 then
        return false,"normalization_map",0,0,0;
    end if;

    E := HyperellipticCurve(g);
    J := Jacobian(E);
    D := J![U,Vn];

    // Locate the reduced elliptic point representing D.  This also catches
    // the case in which the reduction line meets a node preimage.
    Dpt := E![1,0,0];
    foundD := D eq J!0;
    if not foundD then
        for pt in Points(E) do
            if PointClass(J,E,pt) eq D then Dpt:=pt; foundD:=true; break; end if;
        end for;
    end if;
    if not foundD then return false,"elliptic_point_lookup",0,0,0; end if;
    if Dpt[3] ne 0 and Dpt[1]/Dpt[3] eq node then
        return false,"reduced_class_is_node",Order(D),0,0;
    end if;

    eta := SquareRoot(L!Evaluate(g,node));
    vr := L!Evaluate(Vn,node);
    if eta-vr eq 0 or -eta-vr eq 0 then
        return false,"reduction_line_hits_node",Order(D),0,0;
    end if;
    lambdaD := (eta-vr)/(-eta-vr);
    split := IsSquare(Evaluate(g,node));

    ordinary_halves := 0; torus_halves := 0;
    for Qpt in Points(E) do
        H := PointClass(J,E,Qpt);
        if 2*H ne D then continue; end if;
        ordinary_halves +:= 1;

        // Cocycle for adding (Q-infinity) to itself.  For a non-2-torsion
        // point it is the tangent line divided by the vertical through 2Q;
        // the vertical has equal values at the node preimages and cancels
        // from the ratio.  At infinity or at a 2-torsion point the ratio is 1.
        cocycle := L!1;
        if Qpt[3] ne 0 then
            qx := Qpt[1]/Qpt[3];  qy := Qpt[2]/Qpt[3]^2;
            if qy ne 0 then
                slope := Evaluate(Derivative(g),qx)/(2*qy);
                intercept := qy-slope*qx;
                lineplus := eta-L!(slope*node+intercept);
                lineminus := -eta-L!(slope*node+intercept);
                if lineplus eq 0 or lineminus eq 0 then
                    return false,"half_tangent_hits_node",Order(D),
                           ordinary_halves,torus_halves;
                end if;
                cocycle := lineplus/lineminus;
            end if;
        end if;
        target := lambdaD/cocycle;
        // Split torus: F_p^*, order p-1.  Nonsplit torus: norm-one group,
        // order p+1.  These exponent tests are square tests in the relevant
        // torus, not in the ambient F_{p^2}^*.
        torus_square := split select target^((p-1) div 2) eq 1
                              else target^((p+1) div 2) eq 1;
        if torus_square then torus_halves +:= 1; end if;
    end for;
    return torus_halves gt 0,
           split select "complete_split_node" else "complete_nonsplit_node",
           Order(D),ordinary_halves,torus_halves;
end function;

base_double_square:=0; sheet_boundary:=0; collision_bases:=0;
one_node_bases:=0; multi_node_bases:=0; presentations:=0;
chart_fail:=0; complete_tests:=0; complete_divisible:=0;
node_hit:=0; other_unresolved:=0;
split_tests:=0; nonsplit_tests:=0;
elliptic_projection_halves:=0; torus_compatible_halves:=0;
reason_counts:=AssociativeArray();
signature_counts:=AssociativeArray();
signature_reason_counts:=AssociativeArray();
normalization_order_counts:=AssociativeArray();
projection_reason_counts:=AssociativeArray();
projection_order_counts:=AssociativeArray();
projection_tests:=0; projection_divisible:=0; projection_killed:=0;

procedure Bump(~AA,key)
    if IsDefined(AA,key) then AA[key]+:=1; else AA[key]:=1; end if;
end procedure;

print "TARGET_22224_P31_BOUNDARY_STABLE_START", "p",p;

for A in K do if A eq 0 then continue; end if;
for B in K do if B eq 0 then continue; end if;
    C := 1/(A*B);
    R := A^2+B^2+C^2-3;
    S := 1/A^2+1/B^2+1/C^2-3;
    if not IsSquare(R) or not IsSquare(S) then continue; end if;
    base_double_square +:= 1;
    if R eq 0 or S eq 0 then sheet_boundary +:= 1; continue; end if;
    roots2 := [A^2,B^2,C^2,S/R];
    eqs := [<i,j>:i,j in [1..4]|i lt j and roots2[i] eq roots2[j]];
    if #eqs eq 0 then continue; end if;
    collision_bases +:= 1;
    sig := Join([IntegerToString(e[1]) cat "=" cat IntegerToString(e[2]):e in eqs],"+");
    Bump(~signature_counts,sig);
    if #eqs eq 1 then one_node_bases +:= 1; else multi_node_bases +:= 1; end if;

    for rho in RootsBoth(R) do for sigma in RootsBoth(S) do
        ss:=1/(2*rho); tt:=ss^2; uu:=2*tt;
        mags:=[ss*A,ss*B,ss*C]; dd:=2*ss^2*sigma;
        vv:=3*tt^2+dd^2/2;
        triples:=[<mags[1],mags[2],mags[3]>,
                 <mags[1],mags[3],mags[2]>,
                 <mags[2],mags[3],mags[1]>];
        for pairing in [1..3] do
            tr:=triples[pairing];
            ok,f,U,V,reason:=MarkedD12Data(tr[1],tr[2],tr[3],dd,uu,tt,vv);
            presentations +:= 1;
            sing0:=GCD(f,Derivative(f)); repeated0:=Roots(sing0);
            if #repeated0 eq 1 and Degree(sing0) eq 1 then
                pok,pwhy,pord,pnhalf:=ProjectionD3D4(
                    tr[1],tr[2],tr[3],dd,uu,tt,vv,repeated0[1][1]);
                Bump(~projection_reason_counts,pwhy);
                if pok then
                    projection_tests+:=1;
                    Bump(~projection_order_counts,IntegerToString(pord));
                    if pnhalf gt 0 then projection_divisible+:=1;
                    else projection_killed+:=1; end if;
                end if;
            end if;
            if not ok then
                chart_fail+:=1; Bump(~reason_counts,reason);
                Bump(~signature_reason_counts,sig cat "|" cat reason);
                continue;
            end if;
            sing:=GCD(f,Derivative(f));
            repeated:=Roots(sing);
            if #repeated ne 1 or repeated[1][2] ne 1 or Degree(sing) ne 1 then
                other_unresolved+:=1; Bump(~reason_counts,"not_one_node");
                Bump(~signature_reason_counts,sig cat "|not_one_node");
                continue;
            end if;
            node:=repeated[1][1];
            divisible,why,ordD,nhalf,ntorus:=StableDivisible(f,U,V,node);
            Bump(~reason_counts,why);
            Bump(~signature_reason_counts,sig cat "|" cat why);
            if why eq "marked_hits_node" then node_hit+:=1; continue; end if;
            if why[1..Minimum(#why,8)] ne "complete" then
                other_unresolved+:=1; continue;
            end if;
            complete_tests+:=1;
            Bump(~normalization_order_counts,IntegerToString(ordD));
            if why eq "complete_split_node" then split_tests+:=1;
            else nonsplit_tests+:=1; end if;
            elliptic_projection_halves+:=nhalf;
            torus_compatible_halves+:=ntorus;
            if divisible then complete_divisible+:=1; end if;
        end for;
    end for; end for;
end for; end for;

print "P31_BASE_COUNTS",
      "double_square",base_double_square,
      "sheet_boundary",sheet_boundary,
      "collision",collision_bases,
      "one_node",one_node_bases,"multi_node",multi_node_bases;
print "P31_PRESENTATION_COUNTS",
      "presentations",presentations,"chart_fail",chart_fail,
      "complete_generalized_tests",complete_tests,
      "complete_divisible",complete_divisible,
      "node_hit",node_hit,"other_unresolved",other_unresolved;
print "P31_GENERALIZED_DETAIL",
      "split_tests",split_tests,"nonsplit_tests",nonsplit_tests,
      "ordinary_elliptic_halves",elliptic_projection_halves,
      "torus_compatible_halves",torus_compatible_halves;
print "P31_SIGNATURE_COUNTS",Sort([<k,signature_counts[k]>:k in Keys(signature_counts)]);
print "P31_REASON_COUNTS",Sort([<k,reason_counts[k]>:k in Keys(reason_counts)]);
print "P31_SIGNATURE_REASON_COUNTS",
      Sort([<k,signature_reason_counts[k]>:k in Keys(signature_reason_counts)]);
print "P31_NORMALIZATION_ORDER_COUNTS",
      Sort([<k,normalization_order_counts[k]>:k in Keys(normalization_order_counts)]);
print "P31_DIRECT_D3D4_PROJECTION",
      "tests",projection_tests,"ordinary_divisible",projection_divisible,
      "ordinary_killed",projection_killed,
      "orders",Sort([<k,projection_order_counts[k]>:k in Keys(projection_order_counts)]),
      "reasons",Sort([<k,projection_reason_counts[k]>:k in Keys(projection_reason_counts)]);
print "SCOPE complete_generalized_tests are exact Pic0 tests on one-node stable fibres; node_hit, sheet_boundary, and multi-node presentations require blow-up/valuation charts";
print "TARGET_22224_P31_BOUNDARY_STABLE_DONE";
UnsetLogFile();
quit;
