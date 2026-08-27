//////////////////////////////////////////////////////////////////////
// Pairwise elliptic quotient ranks for the seven new repeated-triple
// full-cover fibers found in the B=2000 transverse box search.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);

if not assigned log_file then
    log_file:="results/target_22224_repeated_fibers_quotient_ranks.log";
end if;
SetLogFile(log_file : Overwrite:=true);

fibers := [
 <"F1",[Q!-1470,-630,336],Q!25,Q!5>,
 <"F2",[Q!-720,20,300],Q!-363,Q!17/11>,
 <"F3",[Q!-612,34,289],Q!-338,Q!25/13>,
 <"F4",[Q!-126,28,49],Q!-50,Q!17/5>,
 <"F5",[Q!-112,14,49],Q!-50,Q!13/5>,
 <"F6",[Q!-50,30,45],Q!-48,Q!2>,
 <"F7",[Q!-18,1,16],Q!-50,Q!29/5>
];

function BoundaryProfile(vals,ell)
    F:=GF(ell);
    sq:=[F!0] cat [(F!z)^2:z in vals];
    blocks:=[];
    used:={Integers()|};
    for i in [1..#sq] do
        if i in used then continue; end if;
        block:=[j:j in [i..#sq]|sq[j] eq sq[i]];
        for j in block do Include(~used,j); end for;
        if #block gt 1 then Append(~blocks,block); end if;
    end for;
    return blocks;
end function;

function ExactSquare(x)
    x:=Q!x; if x lt 0 then return false; end if;
    return IsSquare(Z!Numerator(x)) and IsSquare(Z!Denominator(x));
end function;

print "REPEATED_FIBER_RANKS_START","fibers",#fibers;
for fam in fibers do
    name:=fam[1]; fixed:=fam[2]; d0:=fam[3]; t1:=fam[4];
    vals0:=fixed cat [d0]; vals1:=fixed cat [d0*t1^2];
    print "FIBER_START",name,"fixed",fixed,"d0",d0,"t1",t1,
          "member0",vals0,"member1",vals1;
    for tt in [Q!1,t1] do
        vals:=fixed cat [d0*tt^2];
        print "BOUNDARY",name,"t",tt,"p11",BoundaryProfile(vals,11),
              "p13",BoundaryProfile(vals,13);
    end for;
    for ij in [<1,2>,<1,3>,<2,3>] do
        i,j:=Explode(ij); x:=fixed[i]; y:=fixed[j];
        Rx:=(x+d0*T^2)/(x+d0); Ry:=(y+d0*T^2)/(y+d0);
        assert ExactSquare(Evaluate(Rx*Ry,t1));
        try
            C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
            E,phi:=EllipticCurve(C,P0);
            Emin,minmap:=MinimalModel(E);
            lo,hi:=RankBounds(Emin);
            tors:=TorsionSubgroup(Emin);
            print "QUOTIENT",name,ij,"rank_bounds",lo,hi,
                  "torsion",Invariants(tors),"minimal",Emin;
            // Generators are computed only for the most promising fibers in
            // the subsequent MW-sieve job; doing all 21 here is much slower.
        catch err
            print "QUOTIENT_ERROR",name,ij,err`Object;
        end try;
    end for;
    print "FIBER_DONE",name;
end for;
print "REPEATED_FIBER_RANKS_DONE";
UnsetLogFile(); quit;
