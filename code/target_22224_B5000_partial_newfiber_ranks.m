//////////////////////////////////////////////////////////////////////
// Pairwise elliptic quotient ranks for repeated-triple fibres first
// appearing in the live B=5000 transverse full-cover stream.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);

if not assigned log_file then
    log_file:="results/target_22224_B5000_partial_newfiber_ranks.log";
end if;
SetLogFile(log_file : Overwrite:=true);

fibers := [
 <"N1",[Q!-3,10,18],Q!-48/5,Q!2/5>,
 <"N2",[Q!-25,-1,40],Q!4205/98,Q!49/29>,
 <"N3",[Q!-23,-4,50],Q!578/23,Q!37/17>,
 <"N4",[Q!-36,-3,52],Q!768/13,Q!5/2>,
 <"N5",[Q!-31,-2,64],Q!1250/31,Q!41/25>,
 <"N6",[Q!-3,-2,18],Q!588/169,Q!13/7>,
 <"N7",[Q!-15,16,24],Q!-72/5,Q!3/5>,
 <"N8",[Q!-14,-11,35],Q!1375/98,Q!7/5>,
 <"N9",[Q!-46,47,49],Q!-1081/50,Q!5/7>,
 <"N10",[Q!-45,-12,50],Q!3630/49,Q!49/11>,
 // This fibre contains the new local-contact near miss
 // (8,960,-1800,-2535) at T=1/13.
 <"N11",[Q!-64,120,169],Q!-1352/15,Q!1/13>,
 <"N12",[Q!-28,50,55],Q!-3610/77,Q!1/19>,
 // Projectively repeated B=2000 fibres that were invisible to the earlier
 // exact-coordinate grouping.
 <"O1",[Q!-1,4,6],Q!-98/27,Q!3/7>,
 <"O2",[Q!-1,6,9],Q!-507/98,Q!7/13>,
 <"O3",[Q!-3,-2,12],Q!162/49,Q!7/3>,
 <"O4",[Q!-9,10,15],Q!-75/8,Q!1/2>,
 <"O5",[Q!-15,-10,16],Q!50/3,Q!5/3>,
 <"O6",[Q!-15,-1,25],Q!605/27,Q!21/11>,
 <"O7",[Q!-9,-5,30],Q!75/8,Q!5/2>,
 <"O8",[Q!-22,55,70],Q!-1372/25,Q!5/7>,
 <"O9",[Q!1,55,99],Q!5/9,Q!15>,
 // New projective fibres from the live B=10000, a<=1000 stream.
 <"M1",[Q!-16,-5,30],Q!50/3,Q!11/3>,
 <"M2",[Q!-1,25,40],Q!-605/32,Q!17/22>,
 <"M3",[Q!-25,-1,65],Q!405/13,Q!37/27>,
 <"M4",[Q!-49,-2,100],Q!3362/49,Q!61/41>,
 <"M5",[Q!99,244,4026],Q!6,Q!29>
];

function ExactSquare(x)
    x:=Q!x; if x lt 0 then return false; end if;
    return IsSquare(Z!Numerator(x)) and IsSquare(Z!Denominator(x));
end function;

print "B5000_PARTIAL_NEWFIBER_RANKS_START","fibers",#fibers;
for fam in fibers do
    name:=fam[1]; fixed:=fam[2]; d0:=fam[3]; t1:=fam[4];
    print "FIBER_START",name,"fixed",fixed,"d0",d0,"t1",t1,
          "member0",fixed cat [d0],"member1",fixed cat [d0*t1^2];
    for ij in [<1,2>,<1,3>,<2,3>] do
        i,j:=Explode(ij); x:=fixed[i]; y:=fixed[j];
        Rx:=(x+d0*T^2)/(x+d0); Ry:=(y+d0*T^2)/(y+d0);
        assert ExactSquare(Evaluate(Rx,t1));
        assert ExactSquare(Evaluate(Ry,t1));
        try
            C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
            E,phi:=EllipticCurve(C,P0);
            Emin,minmap:=MinimalModel(E);
            lo,hi:=RankBounds(Emin);
            tors:=TorsionSubgroup(Emin);
            print "QUOTIENT",name,ij,"rank_bounds",lo,hi,
                  "torsion",Invariants(tors),"minimal",Emin;
        catch err
            print "QUOTIENT_ERROR",name,ij,err`Object;
        end try;
    end for;
    print "FIBER_DONE",name;
end for;
print "B5000_PARTIAL_NEWFIBER_RANKS_DONE";
UnsetLogFile(); quit;
