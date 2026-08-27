//////////////////////////////////////////////////////////////////////
// Revised target [2,2,4,8]: arithmetic of the three simultaneous-pair
// genus-one fibers through the normalized HLP point.
//
// The fibers force two of F1,F2,F4 on the D-square + F0-square surface.
// This script:
//   * builds the reciprocal binary quartic and its elliptic model;
//   * minimizes the degree-4 model and elliptic Jacobian;
//   * computes torsion, root number, optional conditional-on-GRH RankBounds, and
//     (optionally) a Mordell-Weil group;
//   * maps known/MW points back to the slope model and searches bounded
//     combinations;
//   * imposes C=c^2 and every remaining F_i square condition;
//   * exact-tests J(Q)_tors and applies the n=2,...,12 root-power test.
//
// Typical bounded run (from torsion_jac):
//   magma -b coeff_bound:=12 do_rank:=false do_mw:=false \
//       code/target_2248_hlp_pair_fibers.m
// RankBounds on these large models is class-group dominated, so it is opt-in.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned pair_index then
    pair_index := 0;
elif Type(pair_index) eq MonStgElt then
    pair_index := StringToInteger(pair_index);
end if;

if not assigned coeff_bound then
    coeff_bound := 4;
elif Type(coeff_bound) eq MonStgElt then
    coeff_bound := StringToInteger(coeff_bound);
end if;

if not assigned do_rank then
    do_rank := false;
elif Type(do_rank) eq MonStgElt then
    do_rank := do_rank in {"true","True","1","yes"};
end if;

if not assigned do_mw then
    do_mw := false;
elif Type(do_mw) eq MonStgElt then
    do_mw := do_mw in {"true","True","1","yes"};
end if;

if not assigned use_grh then
    use_grh := true;
elif Type(use_grh) eq MonStgElt then
    use_grh := use_grh in {"true","True","1","yes"};
end if;

if use_grh then
    SetClassGroupBounds("GRH");
end if;

if not assigned write_log then
    write_log := false;
elif Type(write_log) eq MonStgElt then
    write_log := write_log in {"true","True","1","yes"};
end if;
if not assigned log_file then
    log_file := "results/target_2248_hlp_pair_fibers.log";
end if;
if write_log then SetLogFile(log_file : Overwrite := true); end if;

Q := Rationals();
Z := Integers();
R<m> := PolynomialRing(Q);
P<x> := PolynomialRing(Q);

function BoolString(b)
    return b select "true" else "false";
end function;

function IsNonzeroSquareQ(q)
    if q eq 0 then return false,Q!0; end if;
    num := Z!Numerator(q); den := Z!Denominator(q);
    if num lt 0 then return false,Q!0; end if;
    okn,rn := IsSquare(num); okd,rd := IsSquare(den);
    if not okn or not okd then return false,Q!0; end if;
    return true,Q!rn/Q!rd;
end function;

function RatHeight(q)
    return Max(Abs(Z!Numerator(q)),Z!Denominator(q));
end function;

function TargetConicAB(rho,tau,tgt)
    if tgt eq "F1" then
        return (1+rho)*(1+tau),rho/tau;
    elif tgt eq "F2" then
        return (1+rho)*(rho+tau),Q!1/tau;
    end if;
    assert tgt eq "F4";
    return (1+tau)*(rho+tau),rho/tau^2;
end function;

function PairQuartic(A1,B1,A2,B2,q0,y0)
    K := A1*B1;
    den := m^2-K;
    num := q0*m^2-2*y0*m+q0*K;
    return A2*(den^2+B2*num^2);
end function;

function QFromSlope(A,B,q0,y0,s)
    K := A*B;
    den := s^2-K;
    if den eq 0 then return false,Q!0; end if;
    q := (q0*s^2-2*y0*s+q0*K)/den;
    return q ne 0,q;
end function;

function CoverValues(rho,sigma,tau)
    return rho*sigma*tau,
           (1+rho)*(1+sigma)*(1+tau),
           rho*(1+rho)*(rho+sigma)*(rho+tau),
           sigma*(1+sigma)*(rho+sigma)*(sigma+tau),
           tau*(1+tau)*(rho+tau)*(sigma+tau);
end function;

function IntegralSquareTwistModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L,Denominator(Coefficient(f,i)));
    end for;
    return P!(L^2*f);
end function;

function StrictRootPowerCertificate(f)
    C := HyperellipticCurve(f);
    for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
               67,71,73,79,83,89,97,101,103,107,109,113,127,131] do
        try
            fp := ChangeRing(f,GF(pp));
            if Degree(fp) ne Degree(f) or Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C,GF(pp)));
            cs := [Z!Coefficient(Lp,i) : i in [0..4]];
            S<T> := PolynomialRing(Q);
            chi := &+[Q!cs[5-i]*T^i : i in [0..4]];
            if not IsIrreducible(chi) then continue; end if;
            Kf<pi> := NumberField(chi);
            degrees := [Degree(MinimalPolynomial(pi^n)) : n in [2..12]];
            if &and[d eq 4 : d in degrees] then
                return true,pp,chi,degrees;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,0,P!0,[];
end function;

function TorsionPoints(E)
    G,mp := TorsionSubgroup(E);
    return [mp(g) : g in G],Invariants(G);
end function;

function SmallOrder(Pt,B)
    if Pt eq Parent(Pt)!0 then return 1; end if;
    for n in [2..B] do
        if n*Pt eq Parent(Pt)!0 then return n; end if;
    end for;
    return 0;
end function;

StatsFmt := recformat<
    seen_sigma,seen_tuple,pair_points,c_square,full_cover,exact_target,
    simple_hits
>;

procedure TestQ(~stats,label,source,q,rho,tau,dH)
    sigma := q^2*rho/tau;
    skey := Sprint(sigma);
    if skey in stats`seen_sigma then return; end if;
    Include(~stats`seen_sigma,skey);
    stats`pair_points +:= 1;
    if sigma eq 0 or sigma^2 eq 1 or sigma in {-1,-rho,-tau} then return; end if;

    Cval := (sigma^2-rho^2)/(sigma^2-1);
    okC,c := IsNonzeroSquareQ(Cval);
    if not okC then return; end if;
    tup := [rho,Q!1,c,dH];
    if #Set([u^2 : u in tup]) ne 4 then return; end if;
    key := Sprint(Sort([u^2 : u in tup]));
    if key in stats`seen_tuple then return; end if;
    Include(~stats`seen_tuple,key);
    stats`c_square +:= 1;

    F0,F1,F2,F3,F4 := CoverValues(rho,sigma,tau);
    oks := [];
    for F in [F0,F1,F2,F3,F4] do
        okF,rF := IsNonzeroSquareQ(F);
        Append(~oks,okF);
    end for;
    allok := &and oks;
    print "C_SQUARE",label,"source",source,"q_height",RatHeight(q),
          "all_cover",BoolString(allok),"q",q;
    if not allok then return; end if;
    stats`full_cover +:= 1;

    f := x*(x+rho^2)*(x+1)*(x+c^2)*(x+dH^2);
    fI := IntegralSquareTwistModel(f);
    C2 := HyperellipticCurve(fI);
    TG,tmp := TorsionSubgroup(Jacobian(C2));
    invs := Invariants(TG);
    if invs eq [2,2,4,8] then stats`exact_target +:= 1; end if;
    simple,pcert,chi,degrees := StrictRootPowerCertificate(fI);
    if simple then stats`simple_hits +:= 1; end if;
    print "FULL_LIFT",label,"source",source,"tuple",tup;
    print "  CURVE_Y2",fI;
    print "  EXACT_TORSION",invs,"order",#TG;
    print "  ROOT_POWER",BoolString(simple),"p",pcert,"chi",chi,
          "degrees",degrees;
end procedure;

pairs := [<"F1","F2">,<"F1","F4">,<"F2","F4">];

rhoH := Q!58466134224/Q!53109477625;
sigmaH := Q!719363573659505664/Q!749082246897952705;
tauH := Q!307598400/Q!352612321;
cH := Q!58466134224/Q!30294861575;
dH := Q!72946054224/Q!53109477625;
qH := tauH;

print "TARGET_2248_HLP_PAIR_FIBERS";
print "pair_index",pair_index,"coeff_bound",coeff_bound,
      "do_rank",BoolString(do_rank),"do_mw",BoolString(do_mw),
      "use_grh",BoolString(use_grh);
print "HLP",rhoH,sigmaH,tauH,cH,dH;

grand_full := 0;
grand_target := 0;
grand_simple := 0;

for idx in [1..3] do
    if pair_index ne 0 and pair_index ne idx then continue; end if;
    first := pairs[idx][1]; second := pairs[idx][2];
    label := first cat "/" cat second;
    print "PAIR_START",idx,label;

    A1,B1 := TargetConicAB(rhoH,tauH,first);
    A2,B2 := TargetConicAB(rhoH,tauH,second);
    ok1,y1 := IsNonzeroSquareQ(A1*(1+B1*qH^2));
    ok2,y2 := IsNonzeroSquareQ(A2*(1+B2*qH^2));
    assert ok1 and ok2;
    K := A1*B1;
    f4 := PairQuartic(A1,B1,A2,B2,qH,y1);
    tangent := K*qH/y1;
    W0 := y2*(tangent^2-K);
    C4 := HyperellipticCurve(f4);
    Pbase := C4![tangent,W0,1];
    assert Pbase in C4;
    E,phi := EllipticCurve(C4,Pbase);
    Emin := MinimalModel(E);

    print "RECIPROCAL",Coefficient(f4,0) eq Coefficient(f4,4)*K^2,
          Coefficient(f4,1) eq Coefficient(f4,3)*K;
    print "ELLIPTIC_MINIMAL_AINVARIANTS",aInvariants(Emin);
    print "ELLIPTIC_DISCRIMINANT",Discriminant(Emin);
    torsE,tinvs := TorsionPoints(E);
    print "ELLIPTIC_TORSION",tinvs;
    print "ELLIPTIC_ROOT_NUMBER",RootNumber(Emin);

    try
        G1 := GenusOneModel(f4);
        Gmin,trmin,bad := Minimise(G1);
        print "DEGREE4_MINIMISE_BAD_PRIMES",bad;
        print "DEGREE4_MINIMIZED",Eltseq(Gmin);
        try
            Gred,trred := Reduce(Gmin);
            print "DEGREE4_REDUCED",Eltseq(Gred);
        catch er
            print "DEGREE4_REDUCE_ERROR",er`Object;
        end try;
        assert IsIsomorphic(MinimalModel(Jacobian(Gmin)),Emin);
    catch er
        print "DEGREE4_MODEL_ERROR",er`Object;
    end try;

    Pneg := phi(C4![tangent,-W0,1]);
    Pzero := phi(C4![Q!0,-K*y2,1]);
    print "KNOWN_POINT_NEG_ORDER_LE_64",SmallOrder(Pneg,64),"point",Pneg;
    print "KNOWN_POINT_ZERO_ORDER_LE_64",SmallOrder(Pzero,64),"point",Pzero;
    print "KNOWN_POINT_DIFFERENCE_ORDER_LE_64",SmallOrder(Pneg-Pzero,64);
    try
        print "KNOWN_HEIGHT_PAIRING",HeightPairingMatrix([Pneg,Pzero]);
    catch er
        print "KNOWN_HEIGHT_PAIRING_ERROR",er`Object;
    end try;

    if do_rank then
        try
            lo,hi := RankBounds(Emin);
            print "RANK_BOUNDS",lo,hi,"conditional_GRH",use_grh;
        catch er
            print "RANK_BOUNDS_ERROR",er`Object;
        end try;
    end if;

    freepts := [<"Pneg",Pneg>,<"Pzero",Pzero>];
    if do_mw then
        try
            MW,toMin := MordellWeilGroup(Emin);
            print "MORDELL_WEIL_INVARIANTS",Invariants(MW);
            iso := Isomorphism(Emin,E);
            for j in [1..Ngens(MW)] do
                if Order(MW.j) eq 0 then
                    Pj := iso(toMin(MW.j));
                    Append(~freepts,<Sprintf("MW%o",j),Pj>);
                    print "MORDELL_WEIL_GENERATOR",j,Pj;
                end if;
            end for;
        catch er
            print "MORDELL_WEIL_ERROR",er`Object;
        end try;
    end if;

    stats := rec<StatsFmt |
        seen_sigma := {},seen_tuple := {},pair_points := 0,c_square := 0,
        full_cover := 0,exact_target := 0,simple_hits := 0>;
    TestQ(~stats,label,"HLP",qH,rhoH,tauH,dH);

    psi := Inverse(phi);
    P1 := freepts[1][2];
    P2 := freepts[2][2];
    seenE := {};
    for aa in [-coeff_bound..coeff_bound] do
        for bb in [-coeff_bound..coeff_bound] do
            for Tpt in torsE do
                Pe := aa*P1+bb*P2+Tpt;
                ekey := Sprint(Pe);
                if ekey in seenE then continue; end if;
                Include(~seenE,ekey);
                try
                    Pc := psi(Pe);
                    if Pc[3] eq 0 then continue; end if;
                    slope := Q!(Pc[1]/Pc[3]);
                    okq,q := QFromSlope(A1,B1,qH,y1,slope);
                    if not okq then continue; end if;
                    oksecond,tmpy := IsNonzeroSquareQ(A2*(1+B2*q^2));
                    if not oksecond then continue; end if;
                    TestQ(~stats,label,Sprintf("%oPneg+%oPzero+T",aa,bb),
                          q,rhoH,tauH,dH);
                catch er
                    continue;
                end try;
            end for;
        end for;
    end for;

    print "PAIR_DONE",label,"elliptic_points_tested",#seenE,
          "pair_q",stats`pair_points,"c_square",stats`c_square,
          "full_cover",stats`full_cover,"exact_2248",stats`exact_target,
          "simple_hits",stats`simple_hits;
    grand_full +:= stats`full_cover;
    grand_target +:= stats`exact_target;
    grand_simple +:= stats`simple_hits;
end for;

print "SUMMARY", "full_cover",grand_full,"exact_2248",grand_target,
      "simple_hits",grand_simple;
print "TARGET_2248_HLP_PAIR_FIBERS_DONE";
quit;
