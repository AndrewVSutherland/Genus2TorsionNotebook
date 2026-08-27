//////////////////////////////////////////////////////////////////////
// Repeated-support lane on the compact M(12) surface.
//
//   F = L*(L*H^2 + 4*b*(1+x)^2*(w*L-x^2)),
//   L = b+(2*b-1)*x, H=x+w*(1+b*x).
//
// If q=(x-u)^2 and
//
//   A^2-B^2*F=q^5,  deg(A)=5 monic, deg(B)<=2,
//
// with B(u)F(u) nonzero, then v=-A/B mod q is the tangent lift at a
// non-Weierstrass point P and [q,v]=2(P-infinity).  The norm identity
// has divisor 10(P-infinity) on the chosen branch, so P-infinity has
// order dividing 10.  The new lane is exact order 10, equivalently
//
//   5(P-infinity) = W-infinity,
//
// where W is the visible Weierstrass point L=0 (when the rational
// 2-rank is one).
//
// Modes:
//   local       direct finite-field Jacobian census; no norm heuristic;
//   summary     build the exact full quadratic-B norm equations;
//   geometry    bounded modular saturation/decomposition of those equations.
//
// Examples:
//   magma -b mode:="local" p:=7 code/m12_repeated_q_order60.m
//   magma -b mode:="geometry" p:=7 MemGB:=3 \
//       code/m12_repeated_q_order60.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "local"; end if;
require mode in {"local","summary","geometry"}:
        "mode must be local, summary, or geometry";
if not assigned p then p := 7;
elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned MemGB then MemGB := 3;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned sample_limit then sample_limit := 12;
elif Type(sample_limit) eq MonStgElt then
    sample_limit := StringToInteger(sample_limit);
end if;
if not assigned do_primary then do_primary := true;
elif Type(do_primary) eq MonStgElt then
    do_primary := do_primary in {"true","True","1"};
end if;
SetMemoryLimit(MemGB*10^9);

function ForcedA(F,q,B)
    P := Parent(F); x := P.1;
    A := x^5;
    for d in [9,8,7,6,5] do
        E := A^2-B^2*F-q^5;
        A +:= (-Coefficient(E,d)/2)*x^(d-5);
    end for;
    E := A^2-B^2*F-q^5;
    assert &and[Coefficient(E,d) eq 0:d in [5..10]];
    return A,E;
end function;

function NormData(k : order := "grevlex")
    R<b,w,u,b0,b1,b2> := PolynomialRing(k,6,order);
    K := FieldOfFractions(R);
    P<x> := PolynomialRing(K);
    L := b+(2*b-1)*x;
    H := x+w*(1+b*x);
    F := L*(L*H^2+4*b*(1+x)^2*(w*L-x^2));
    q := (x-u)^2;
    B := b0+b1*x+b2*x^2;
    A,E := ForcedA(F,q,B);
    eqs := [R!Coefficient(E,i):i in [0..4]];
    assert E eq &+[K!eqs[i+1]*x^i:i in [0..4]];
    Fu := R!Evaluate(F,u);
    Bu := b0+b1*u+b2*u^2;
    discF := R!Discriminant(F);
    cheap_names := ["b","w","b-1","2b-1","Fu","Bu"];
    cheap := [b,w,b-1,2*b-1,Fu,Bu];
    return R,P,F,q,B,A,eqs,cheap_names,cheap,discF;
end function;

function IsUnitIdeal(I)
    return &or[g ne 0 and TotalDegree(g) eq 0:g in Basis(I)];
end function;

procedure Report(label,I)
    B := Basis(I);
    print "IDEAL",label,"unit",IsUnitIdeal(I),"basis_len",#B,
          "degrees",[TotalDegree(f):f in B];
    if IsUnitIdeal(I) then return; end if;
    try
        d,ds := Dimension(I);
        print "IDEAL",label,"dimension",d,"component_degrees",ds;
    catch err
        print "IDEAL",label,"dimension_failed",err`Object;
    end try;
    try print "IDEAL",label,"degree",Degree(I);
    catch err print "IDEAL",label,"degree_failed",err`Object; end try;
end procedure;

if mode eq "local" then
    require IsPrime(p) and p notin {2,3,5}: "use a prime away from 2,3,5";
    k := GF(p); P<x> := PolynomialRing(k);
    base_smooth := 0; point_pairs := 0; exact5 := 0; exact10 := 0;
    relation10 := 0; rank1_2 := 0; samples := [];
    for b in k do for w in k do
        if b*w*(b-1)*(2*b-1) eq 0 then continue; end if;
        L := b+(2*b-1)*x;
        H := x+w*(1+b*x);
        F := L*(L*H^2+4*b*(1+x)^2*(w*L-x^2));
        if Degree(F) ne 5 or Discriminant(F) eq 0 then continue; end if;
        base_smooth +:= 1;
        C := HyperellipticCurve(F); J := Jacobian(C); O := J!0;
        wr := -b/(2*b-1);
        T := J![x-wr,k!0];
        assert T ne O and 2*T eq O;
        fac := Factorization(F div L);
        two_rank_one := (#fac eq 1 and Degree(fac[1][1]) eq 4);
        if two_rank_one then rank1_2 +:= 1; end if;
        for u in k do
            fu := Evaluate(F,u);
            if fu eq 0 or not IsSquare(fu) then continue; end if;
            for y in [Sqrt(fu),-Sqrt(fu)] do
                // In odd characteristic y and -y are distinct since fu != 0.
                E := J![x-u,y];
                point_pairs +:= 1;
                if 5*E eq O then
                    exact5 +:= 1;
                    if #samples lt sample_limit then
                        Append(~samples,<"ord5",b,w,u,y,two_rank_one>);
                    end if;
                elif 5*E eq T then
                    relation10 +:= 1;
                    if 2*E ne O and 10*E eq O then
                        exact10 +:= 1;
                        if #samples lt sample_limit then
                            Append(~samples,<"ord10",b,w,u,y,two_rank_one>);
                        end if;
                    end if;
                end if;
            end for;
        end for;
    end for; end for;
    print "M12_REPEATED_Q_LOCAL","p",p,"smooth_bases",base_smooth,
          "two_rank_one_bases",rank1_2,"nonbranch_point_pairs",point_pairs,
          "exact_order5_points",exact5,"fiveP_equals_T",relation10,
          "exact_order10_points",exact10,"samples",samples;
    quit;
end if;

k := mode eq "summary" select Rationals() else GF(p);
if mode eq "geometry" then
    require IsPrime(p) and p notin {2,3,5}: "use a prime away from 2,3,5";
end if;
R,P,F,q,B,A,eqs,names,factors,discF := NormData(k);
print "M12_REPEATED_Q_NORM_FORMULAS_PASS";
print "EQUATION_SHAPES",[<TotalDegree(f),#Terms(f)>:f in eqs];
print "BOUNDARY_SHAPES",
      [<names[i],TotalDegree(factors[i]),#Terms(factors[i])>:
       i in [1..#factors]],
      <"discF",TotalDegree(discF),#Terms(discF)>;
if mode eq "summary" then quit; end if;

I := ideal<R|eqs>;
Report("raw",I);
J := I;
for i in [1..#factors] do
    nf := NormalForm(factors[i],Basis(J));
    print "SAT_BEGIN",names[i],"degree",TotalDegree(nf),"terms",#Terms(nf);
    time J := Saturation(J,ideal<R|nf>);
    Report(Sprintf("sat_%o",names[i]),J);
    if IsUnitIdeal(J) then break; end if;
end for;
if not IsUnitIdeal(J) then
    nf := NormalForm(discF,Basis(J));
    print "SAT_BEGIN discF","degree",TotalDegree(nf),"terms",#Terms(nf);
    time J := Saturation(J,ideal<R|nf>);
    Report("full_open",J);
end if;
if do_primary and not IsUnitIdeal(J) then
    print "PRIMARY_BEGIN";
    try
        time comps := PrimaryDecomposition(J);
        print "PRIMARY_COUNT",#comps;
        for i in [1..#comps] do
            Report(Sprintf("component_%o",i),comps[i]);
            BC := Basis(comps[i]);
            if #BC le 15 and &and[#Terms(g) le 100:g in BC] then
                print "COMPONENT_BASIS",i,BC;
            end if;
        end for;
    catch err
        print "PRIMARY_FAILED",err`Object;
    end try;
end if;
print "M12_REPEATED_Q_GEOMETRY_DONE";
quit;
