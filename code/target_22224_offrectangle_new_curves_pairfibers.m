//////////////////////////////////////////////////////////////////////
// Seven NEW elliptic pair-scaling fibers mined from the 56 primitive
// off-rectangle tor2228 points.
//
// On a fiber
//
//     (a,b,c,d) = (lambda*a0,b0,c0,lambda*d0),
//
// abcd stays square automatically.  Relative to lambda=1, put
//
//   Ra=((a*l+b)(a*l+c))/((a+b)(a+c)),
//   Rb=((a*l+b)(d*l+b))/((a+b)(d+b)),
//   Rc=((a*l+c)(d*l+c))/((a+c)(d+c)).
//
// The full A2228 conditions are Ra,Rb,Rc squares.  With p=ya/yb and
// q=yc/ya, eliminating lambda gives A(p^2)q^2+B(p^2)=0 and hence the
// genus-one quotient z^2=-A(p^2)B(p^2).  We enumerate its Mordell--Weil
// group and retain exactly those points for which Ra is itself a square;
// these are precisely the rational points lifting to the full fiber.
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
Q:=Rationals(); Z:=Integers(); PR<p>:=PolynomialRing(Q);

if not assigned LatticeBound then LatticeBound:=8; end if;
if Type(LatticeBound) eq MonStgElt then LatticeBound:=StringToInteger(LatticeBound); end if;
if not assigned MultipleBound then MultipleBound:=100; end if;
if Type(MultipleBound) eq MonStgElt then MultipleBound:=StringToInteger(MultipleBound); end if;
if not assigned FiberStart then FiberStart:=1; end if;
if Type(FiberStart) eq MonStgElt then FiberStart:=StringToInteger(FiberStart); end if;
if not assigned FiberEnd then FiberEnd:=7; end if;
if Type(FiberEnd) eq MonStgElt then FiberEnd:=StringToInteger(FiberEnd); end if;
if not assigned PrintPoints then PrintPoints:=false; end if;
if not assigned DoFullMW then DoFullMW:=true; end if;
if Type(DoFullMW) eq MonStgElt then DoFullMW:=DoFullMW eq "true"; end if;
if Type(PrintPoints) eq MonStgElt then PrintPoints:=PrintPoints eq "true"; end if;
if not assigned HeightCeiling then HeightCeiling:=10^40; end if;
if Type(HeightCeiling) eq MonStgElt then HeightCeiling:=StringToInteger(HeightCeiling); end if;
if not assigned output_file then
    output_file:="results/target_22224_offrectangle_new_curves_pairfiber_candidates.tsv";
end if;
if not assigned log_file then
    log_file:="results/target_22224_offrectangle_new_curves_pairfibers.log";
end if;

SetLogFile(log_file:Overwrite:=true);
out:=Open(output_file,"w");
fprintf out,"fiber\tsource\tlambda\ta\tb\tc\td\tprimitive_a\tprimitive_b\tprimitive_c\tprimitive_d\tknown_seed\tthree_survivor\tgcd_bound\n";

prime_list := [11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
               79,83,89,97,101,103,107,109,113,127,131,137,139,149,
               151,157,163,167,173,179,181,191,193,197,199];

// Entries are <name,[moving a,fixed b,fixed c,moving d], second lambda>.
fibers := [
    <"P1",[Q!1,Q!55,Q!99,Q!125],Q!1089/25>,
    <"P2",[Q!20,Q!225,Q!304,Q!380],Q!9>,
    <"P3",[Q!13,Q!276,Q!624,Q!828],Q!16>,
    <"P4",[Q!41,Q!256,Q!800,Q!1312],Q!6400/1681>,
    <"P5",[Q!2,Q!46,Q!292,Q!1679],Q!4>,
    <"P6",[Q!18,Q!158,Q!711,Q!1764],Q!6241/1764>,
    <"P7",[Q!50,Q!791,Q!2800,Q!5650],Q!196/25>
];

function ExactSquare(x)
    x:=Q!x;
    if x lt 0 then return false,Q!0; end if;
    okN,rtN:=IsSquare(Z!Numerator(x)); if not okN then return false,Q!0; end if;
    okD,rtD:=IsSquare(Z!Denominator(x)); if not okD then return false,Q!0; end if;
    return true,Q!rtN/Q!rtD;
end function;

function PrimitiveTuple(vals)
    den:=LCM([Denominator(Q!z):z in vals]);
    ints:=[Z!(Q!z*den):z in vals]; g:=GCD(ints);
    if g ne 0 then ints:=[z div g:z in ints]; end if;
    return Sort([Abs(z):z in ints]);
end function;

function Rectangle(vals)
    a,b,c,d:=Explode(vals);
    return a*b eq c*d or a*c eq b*d or a*d eq b*c;
end function;

function SmoothTuple(vals)
    return not &or[z eq 0:z in vals] and #Set([z^2:z in vals]) eq 4;
end function;

function FullCover(vals)
    a,b,c,d:=Explode(vals);
    rad:=[a*b*c*d,a*(a+b)*(a+c)*(a+d),
         b*(a+b)*(b+c)*(b+d),c*(a+c)*(b+c)*(c+d)];
    return &and[ExactSquare(z):z in rad];
end function;

function ThreeBound(vals)
    f:=p*&*[p+z^2:z in vals]; g:=0; used:=[];
    for ell in prime_list do
        try
            fp:=ChangeRing(f,GF(ell));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            n:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g:=(g eq 0) select n else GCD(g,n); Append(~used,<ell,n,g>);
            if g mod 3 ne 0 then return false,g,used; end if;
        catch e continue;
        end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

function LatticePoints(gens,bound,E)
    pts:={E!0};
    for gen in gens do
        nxt:={};
        for ep in pts do for n in [-bound..bound] do Include(~nxt,ep+n*gen); end for; end for;
        pts:=nxt;
    end for;
    return pts;
end function;

knownSeeds:={
 <1,55,99,125>,<99,125,225,12375>,<20,225,304,380>,<180,225,304,3420>,
 <13,276,624,828>,<52,69,156,3312>,<41,256,800,1312>,<200,328,1025,6400>,
 <2,46,292,1679>,<4,23,146,3358>,<18,158,711,1764>,<79,196,882,7742>,
 <50,791,2800,5650>,<56,113,400,6328>
};

seen:={}; quotients:=0; rankPositive:=0; searched:=0; liftPoints:=0;
newPoints:=0; threeSurvivors:=0; heightSkipped:=0;
print "TARGET_22224_OFFRECTANGLE_PAIRFIBERS_START",
      "LatticeBound",LatticeBound,"MultipleBound",MultipleBound,
      "FiberStart",FiberStart,"FiberEnd",FiberEnd,
      "DoFullMW",DoFullMW,"HeightCeiling",HeightCeiling;

for fi in [Max(1,FiberStart)..Min(#fibers,FiberEnd)] do
    fib:=fibers[fi];
    name:=fib[1]; vals0:=fib[2]; lam1:=fib[3];
    a,b,c,d:=Explode(vals0);
    Kp:=(d+b)/(a+c); Kq:=(a+b)/(d+c);
    Ap := a*(Kp*c-b*p^2)+b*(d*p^2-Kp*a);
    Bp := -Kq*d*(Kp*c-b*p^2)-Kq*c*(d*p^2-Kp*a);
    quartic := -Ap*Bp;
    assert Degree(quartic) eq 4 and Discriminant(quartic) ne 0;
    C:=HyperellipticCurve(quartic);
    Aval0:=Evaluate(Ap,Q!1); Bval0:=Evaluate(Bp,Q!1);
    assert Aval0+Bval0 eq 0;
    P0:=C![Q!1,Aval0,Q!1]; assert P0 in C;

    function Ratios(lam)
        Ra:=((a*lam+b)*(a*lam+c))/((a+b)*(a+c));
        Rb:=((a*lam+b)*(d*lam+b))/((a+b)*(d+b));
        Rc:=((a*lam+c)*(d*lam+c))/((a+c)*(d+c));
        return Ra,Rb,Rc;
    end function;
    Ra1,Rb1,Rc1:=Ratios(lam1);
    oka,ya:=ExactSquare(Ra1); okb,yb:=ExactSquare(Rb1); okc,yc:=ExactSquare(Rc1);
    assert oka and okb and okc;
    p1:=ya/yb; q1:=yc/ya; A1:=Evaluate(Ap,p1);
    assert q1^2*A1+Evaluate(Bp,p1) eq 0;
    P1:=C![p1,q1*A1,Q!1]; assert P1 in C;

    E,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi); Q1:=phi(P1);
    ord:=Order(Q1); if ord eq 0 then rankPositive+:=1; end if;
    Emin:=MinimalModel(E); lo:=-1; hi:=-1;
    TG,tmap:=TorsionSubgroup(E); torsPts:={tmap(g):g in TG};
    invMW:=[]; basis:=[Q1];
    if DoFullMW then
        try lo,hi:=RankBounds(Emin); catch e lo:=-1;hi:=-1; end try;
        MW,mwmap:=MordellWeilGroup(E); invMW:=Invariants(MW);
        basis:=[mwmap(MW.i):i in [1..#invMW] | invMW[i] eq 0];
    end if;
    quotients+:=1;
    print "PAIRFIBER",name,"base",vals0,"lambda1",lam1,
          "quartic",quartic,"known_order",ord,"rank_bounds",lo,hi,
          "MW",invMW,"torsion",Invariants(TG),"minimal",Emin;

    epoints:=LatticePoints(basis,LatticeBound,E);
    for n in [-MultipleBound..MultipleBound] do Include(~epoints,n*Q1); end for;
    withTors:={}; for ep in epoints do for tor in torsPts do Include(~withTors,ep+tor); end for; end for;
    searched+:=#withTors;
    for ep in withTors do
        try cp:=Einv(ep); catch e continue; end try;
        if cp[3] eq 0 then continue; end if;
        pp:=Q!(cp[1]/cp[3]); zz:=Q!(cp[2]/cp[3]^2);
        Aval:=Evaluate(Ap,pp); if Aval eq 0 then continue; end if;
        qq:=zz/Aval;
        den:=pp^2*d-Kp*a; if den eq 0 then continue; end if;
        lam:=(Kp*c-pp^2*b)/den;
        if pp^2 ne (Kp*(a*lam+c)/(d*lam+b)) then continue; end if;
        if qq^2 ne (Kq*(d*lam+c)/(a*lam+b)) then continue; end if;
        Ra,Rb,Rc:=Ratios(lam);
        ok,yaLift:=ExactSquare(Ra); if not ok then continue; end if;
        // The ratio equations now force Rb and Rc to be squares too.
        okb2,ybLift:=ExactSquare(Rb); okc2,ycLift:=ExactSquare(Rc);
        if not okb2 or not okc2 then continue; end if;
        vals:=[a*lam,b,c,d*lam];
        if not SmoothTuple(vals) or not FullCover(vals) then continue; end if;
        prim:=PrimitiveTuple(vals); key:=<prim[1],prim[2],prim[3],prim[4]>;
        if key in seen then continue; end if; Include(~seen,key); liftPoints+:=1;
        if Maximum(prim) gt HeightCeiling then heightSkipped+:=1; continue; end if;
        known:=key in knownSeeds;
        if not known and not Rectangle(prim) then newPoints+:=1; end if;
        survives,g,used:=ThreeBound(vals); if survives then threeSurvivors+:=1; end if;
        fprintf out,"%o\tMW\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                name,lam,vals[1],vals[2],vals[3],vals[4],
                prim[1],prim[2],prim[3],prim[4],known select 1 else 0,
                survives select 1 else 0,g;
        if PrintPoints or known or survives then
            print "PAIRFIBER_LIFT",name,"lambda",lam,"primitive",prim,
                  "known",known,"rectangle",Rectangle(prim),
                  "three_survivor",survives,"gcd",g,"used",used;
        else
            print "PAIRFIBER_LIFT_SUMMARY",name,
                  "height_digits",#IntegerToString(Maximum(prim)),
                  "three_survivor",survives,"gcd",g;
        end if;
        if survives then
            f:=p*&*[p+z^2:z in vals]; G,jmap:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
            print "EXACT_TORSION",name,prim,Invariants(G);
        end if;
    end for;
end for;

delete out;
print "TARGET_22224_OFFRECTANGLE_PAIRFIBERS_DONE",
      "quotients",quotients,"known_points_infinite",rankPositive,
      "searched_E_points",searched,"full_lifts",liftPoints,
      "height_skipped",heightSkipped,"new_offrectangle_full",newPoints,
      "three_survivors",threeSurvivors;
UnsetLogFile(); quit;
