//////////////////////////////////////////////////////////////////////
// Elliptic-quotient mining on the three repeated-coordinate fibers in
// the primitive off-rectangle tor2228 bank.
//
// Fix (a,b,c) and write d=d0*T^2.  Relative to the known point T=1,
// the full cover is the genus-5 multiquadratic curve
//
//   Y_a^2=(a+d0*T^2)/(a+d0), etc.
//
// Every pair gives an elliptic quotient
//
//   Z^2=R_i(T)R_j(T).
//
// The second bank point gives T=T1 on every quotient.  We compute ranks,
// enumerate multiples and torsion translates, then impose all three
// individual square conditions again.  New full-cover tuples are subjected
// to the safe direct point-count 3-primary filter.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);
Q := Rationals(); Z := Integers();
R<T> := PolynomialRing(Q);

if not assigned MultipleBound then MultipleBound := 100; end if;
if Type(MultipleBound) eq MonStgElt then MultipleBound := StringToInteger(MultipleBound); end if;
if not assigned PointBound then PointBound := 1000; end if;
if Type(PointBound) eq MonStgElt then PointBound := StringToInteger(PointBound); end if;
if not assigned LatticeBound then LatticeBound := 4; end if;
if Type(LatticeBound) eq MonStgElt then LatticeBound := StringToInteger(LatticeBound); end if;
if not assigned FiberStart then FiberStart := 1; end if;
if Type(FiberStart) eq MonStgElt then FiberStart := StringToInteger(FiberStart); end if;
if not assigned FiberEnd then FiberEnd := 4; end if;
if Type(FiberEnd) eq MonStgElt then FiberEnd := StringToInteger(FiberEnd); end if;
if not assigned output_file then
    output_file := "results/target_22224_offrectangle_new_curves_fiber_candidates.tsv";
end if;
if not assigned log_file then
    log_file := "results/target_22224_offrectangle_new_curves_fibers.log";
end if;

SetLogFile(log_file : Overwrite := true);
out := Open(output_file,"w");
fprintf out,"fiber\tpair\tsource\tt\ta\tb\tc\td\tprimitive_a\tprimitive_b\tprimitive_c\tprimitive_d\tin_bank\tthree_survivor\tgcd_bound\n";

prime_list := [11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
               79,83,89,97,101,103,107,109,113,127,131,137,139,149,
               151,157,163,167,173,179,181,191,193,197,199];

fibers := [
    <"F1",[Q!46,Q!292,Q!1679],Q!2,Q!19>,
    <"F2",[Q!71,Q!3650,Q!10366],Q!2116,Q!37/23>,
    <"F3",[Q!276,Q!624,Q!828],Q!13,Q!11>,
    // Signed transverse-box discovery through the best 3-contact near miss:
    // (50,528,-726,-891) -> (722,528,-726,-891), 722/50=(19/5)^2.
    <"P8",[Q!528,Q!-726,Q!-891],Q!50,Q!19/5>
];

function ExactSquare(x)
    x := Q!x;
    if x lt 0 then return false,Q!0; end if;
    okN,rtN := IsSquare(Z!Numerator(x));
    if not okN then return false,Q!0; end if;
    okD,rtD := IsSquare(Z!Denominator(x));
    if not okD then return false,Q!0; end if;
    return true,Q!rtN/Q!rtD;
end function;

function PrimitiveTuple(vals)
    den := LCM([Denominator(Q!z):z in vals]);
    ints := [Z!(Q!z*den):z in vals];
    g := GCD(ints);
    if g ne 0 then ints := [z div g:z in ints]; end if;
    if &and[z le 0:z in ints] then ints := [-z:z in ints]; end if;
    return Sort([Abs(z):z in ints]);
end function;

function Rectangle(vals)
    a,b,c,d := Explode(vals);
    return a*b eq c*d or a*c eq b*d or a*d eq b*c;
end function;

function SmoothTuple(vals)
    if &or[z eq 0:z in vals] then return false; end if;
    return #Set([z^2:z in vals]) eq 4;
end function;

function FullCover(vals)
    a,b,c,d := Explode(vals);
    rad := [
        a*b*c*d,
        a*(a+b)*(a+c)*(a+d),
        b*(a+b)*(b+c)*(b+d),
        c*(a+c)*(b+c)*(c+d)
    ];
    return &and[ExactSquare(z):z in rad];
end function;

function ThreeBound(vals)
    f := T*&*[T+z^2:z in vals];
    g := 0; used := [];
    for p in prime_list do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            n := Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g := (g eq 0) select n else GCD(g,n);
            Append(~used,<p,n,g>);
            if g mod 3 ne 0 then return false,g,used; end if;
        catch e
            continue;
        end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

function LatticePoints(gens,bound,E)
    pts := {E!0};
    for gen in gens do
        nxt := {};
        for ep in pts do
            for n in [-bound..bound] do Include(~nxt,ep+n*gen); end for;
        end for;
        pts := nxt;
    end for;
    return pts;
end function;

// The 56 primitive off-rectangle tuples, for exact novelty checks.
bank56 := {
    <1,55,99,125>,<20,225,304,380>,<95,125,169,931>,<13,276,624,828>,
    <41,256,800,1312>,<67,80,245,1675>,<2,46,292,1679>,<276,624,828,1573>,
    <46,292,722,1679>,<18,158,711,1764>,<31,558,900,2178>,<95,121,529,2375>,
    <256,959,1096,2744>,<52,69,156,3312>,<4,23,146,3358>,<180,225,304,3420>,
    <217,1008,1519,3481>,<169,1271,2009,3751>,<41,1808,2809,4633>,
    <121,1919,3211,4949>,<833,2529,4496,4913>,<99,244,4026,5046>,
    <17,4208,4471,5329>,<50,791,2800,5650>,<56,113,400,6328>,
    <200,328,1025,6400>,<64,136,225,6664>,<121,423,729,6768>,
    <23,544,4352,6624>,<25,264,936,6864>,<4,3131,6076,6464>,
    <50,2116,3431,6862>,<18,1351,6174,6948>,<79,196,882,7742>,
    <36,439,3042,7902>,<1136,2704,6976,7739>,<21,371,3021,8379>,
    <496,562,2738,8711>,<18,686,4932,8631>,<801,1600,2225,9025>,
    <645,1155,3773,9675>,<1369,1711,2349,9971>,<2,448,2914,10199>,
    <71,2116,3650,10366>,<71,3650,5476,10366>,<26,2140,5350,10985>,
    <99,125,225,12375>,<6,219,3358,13248>,<1271,2480,10000,13120>,
    <333,2450,11328,13098>,<1625,3100,7750,13754>,<529,5191,8381,14499>,
    <46,1682,2692,15479>,<484,5795,9500,15616>,<68,882,1854,15759>,
    <306,5476,11575,15742>
};

seenTuples := {}; totalQuotients:=0; positiveRank:=0; searchedPoints:=0;
fullPoints:=0; newFull:=0; threeSurvivors:=0;

print "TARGET_22224_OFFRECTANGLE_FIBERS_START",
      "MultipleBound",MultipleBound,"PointBound",PointBound,
      "LatticeBound",LatticeBound;

for fibno in [Max(1,FiberStart)..Min(#fibers,FiberEnd)] do
    fib := fibers[fibno];
    name := fib[1]; fixed := fib[2]; d0 := fib[3]; t1 := fib[4];
    print "FIBER_START",name,"fixed",fixed,"d0",d0,"known_t",t1;
    for ij in [<1,2>,<1,3>,<2,3>] do
        i,j := Explode(ij); x:=fixed[i]; y:=fixed[j];
        Rx := (x+d0*T^2)/(x+d0);
        Ry := (y+d0*T^2)/(y+d0);
        quartic := Rx*Ry;
        C := HyperellipticCurve(quartic);
        P0 := C![Q!1,Q!1,Q!1]; assert P0 in C;
        okx,yx := ExactSquare(Evaluate(Rx,t1));
        oky,yy := ExactSquare(Evaluate(Ry,t1));
        assert okx and oky;
        P1 := C![t1,yx*yy,Q!1]; assert P1 in C;
        E,phi := EllipticCurve(C,P0);
        Einv := Inverse(phi);
        Q1 := phi(P1);
        totalQuotients +:= 1;
        ord := Order(Q1);
        if ord eq 0 then positiveRank +:= 1; end if;
        Emin := MinimalModel(E);
        try lo,hi := RankBounds(Emin); catch e lo:=-1; hi:=-1; end try;
        TG,tmap := TorsionSubgroup(E);
        torsPts := {tmap(g):g in TG};
        print "QUOTIENT",name,ij,"quartic",quartic,
              "known_point_order",ord,"rank_bounds",lo,hi,
              "torsion",Invariants(TG),"minimal",Emin;

        generators := {Q1};
        try
            for ep in PointSearch(E,PointBound) do
                if Order(ep) eq 0 then Include(~generators,ep); end if;
            end for;
        catch e
            print "POINTSEARCH_ERROR",name,ij,e`Object;
        end try;
        basis := [];
        try
            MW,mwmap := MordellWeilGroup(E);
            invMW := Invariants(MW);
            basis := [mwmap(MW.i):i in [1..#invMW] | invMW[i] eq 0];
            print "MW_GROUP",name,ij,"invariants",invMW,"basis_size",#basis;
        catch e
            basis := Setseq(generators);
            print "MW_GROUP_ERROR",name,ij,e`Object;
        end try;
        if #basis eq 0 then basis := Setseq(generators); end if;
        print "QUOTIENT_GENERATORS",name,ij,"pointsearch",#generators,
              "lattice_basis",#basis;

        epoints := LatticePoints(basis,LatticeBound,E);
        // Retain a deeper one-generator walk along every found point as well.
        for gen in generators do
            for n in [-MultipleBound..MultipleBound] do Include(~epoints,n*gen); end for;
        end for;
        withTors := {};
        for ep in epoints do for tor in torsPts do Include(~withTors,ep+tor); end for; end for;
        epoints := withTors;
        searchedPoints +:= #epoints;
        for ep in epoints do
            try cp := Einv(ep); catch e continue; end try;
            if cp[3] eq 0 then continue; end if;
            tt := Q!(cp[1]/cp[3]);
            dd := d0*tt^2;
            vals := [fixed[1],fixed[2],fixed[3],dd];
            if not SmoothTuple(vals) then continue; end if;
            allsq := true;
            for z in fixed do
                ok,sq := ExactSquare((z+dd)/(z+d0));
                if not ok then allsq:=false; break; end if;
            end for;
            if not allsq or not FullCover(vals) then continue; end if;
            prim := PrimitiveTuple(vals);
            key := <prim[1],prim[2],prim[3],prim[4]>;
            if key in seenTuples then continue; end if;
            Include(~seenTuples,key); fullPoints +:= 1;
            inbank := key in bank56;
            if not inbank and not Rectangle(prim) then newFull +:= 1; end if;
            survives,g,used := ThreeBound(vals);
            if survives then threeSurvivors +:=1; end if;
            fprintf out,"%o\t%o-%o\tMW\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                    name,i,j,tt,vals[1],vals[2],vals[3],vals[4],
                    prim[1],prim[2],prim[3],prim[4],
                    inbank select 1 else 0,survives select 1 else 0,g;
            print "FULL_FIBER_POINT",name,ij,"t",tt,"primitive",prim,
                  "in_bank",inbank,"rectangle",Rectangle(prim),
                  "three_survivor",survives,"gcd",g,"used",used;
            if survives then
                f := T*&*[T+z^2:z in vals];
                G,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
                print "EXACT_TORSION",Invariants(G);
            end if;
        end for;
    end for;
    print "FIBER_DONE",name;
end for;

delete out;
print "TARGET_22224_OFFRECTANGLE_FIBERS_DONE",
      "quotients",totalQuotients,"positive_rank_known_point",positiveRank,
      "searched_E_points",searchedPoints,"full_points",fullPoints,
      "new_offrectangle_full",newFull,"three_survivors",threeSurvivors;
UnsetLogFile();
quit;
