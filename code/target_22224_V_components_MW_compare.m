//////////////////////////////////////////////////////////////////////
// Exact MW-to-T comparison for the surviving V1/V2 rank-two sieve rows.
// This does not rerun the sieve: it maps the recorded coefficient rows
// back to P^1(T), reports local valuations, and checks the three remaining
// full-cover square conditions exactly.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Q := Rationals(); Z := Integers(); R<T> := PolynomialRing(Q);
SetLogFile("results/target_22224_V_components_MW_compare.log":Overwrite:=true);

function Vp(q,p)
    if q eq 0 then return 1000000; end if;
    n:=Numerator(Q!q); d:=Denominator(Q!q); v:=0;
    while n mod p eq 0 do n div:=p; v+:=1; end while;
    while d mod p eq 0 do d div:=p; v-:=1; end while;
    return v;
end function;

procedure Analyze(name,fixed,d0,known_t,pair,rows,primes)
    Rx := (fixed[pair[1]]+d0*T^2)/(fixed[pair[1]]+d0);
    Ry := (fixed[pair[2]]+d0*T^2)/(fixed[pair[2]]+d0);
    C := HyperellipticCurve(Rx*Ry); P0 := C![Q!1,Q!1,Q!1];
    Eraw,phi := EllipticCurve(C,P0); Einv := Inverse(phi);
    E,minmap := MinimalModel(Eraw); minmapinv := Inverse(minmap);
    free := [g:g in Generators(E)|Order(g) eq 0];
    TG,tmap := TorsionSubgroup(E); tors := [tmap(g):g in TG];
    print "MW_COMPARE_START",name,"E",E,"rank",#free,"torsion",Invariants(TG);
    for row in rows do
        m:=row[1]; n:=row[2]; ti:=row[3];
        ep:=m*free[1]+n*free[2]+tors[ti];
        cp:=Einv(minmapinv(ep));
        if cp[3] eq 0 then
            print "MW_ROW",name,row,"T=geometric_infinity";
            continue;
        end if;
        tt:=Q!(cp[1]/cp[3]); dd:=d0*tt^2;
        sq:=[]; vals:=[];
        for z in fixed do
            q:=(z+dd)/(z+d0); Append(~vals,q); Append(~sq,IsSquare(q));
        end for;
        loc:=[];
        for p in primes do
            if Vp(tt,p) ge 0 and Denominator(tt) mod p ne 0 then
                Append(~loc,<p,"finite",Z!(Numerator(tt)*
                    (Z!InverseMod(Denominator(tt) mod p,p)) mod p),Vp(tt,p)>);
            else
                Append(~loc,<p,"infinity",0,Vp(tt,p)>);
            end if;
        end for;
        h:=Max(Abs(Numerator(tt)),Denominator(tt));
        if h lt 10^30 then
            print "MW_ROW",name,row,"T",tt,"height",h,
                  "local",loc,"cover_squares",sq;
        else
            print "MW_ROW",name,row,"T_omitted_digits",
                  <#IntegerToString(Abs(Numerator(tt))),#IntegerToString(Denominator(tt))>,
                  "local",loc,"cover_squares",sq;
        end if;
    end for;
end procedure;

Analyze("V1",[Q!-18,20,75],Q!-1470/121,Q!11/49,<1,2>,
        [<-2,-1,1>],[13,11]);
// All six recorded V2 rows are congruent to (m,n,ti)=(3,3,1) modulo
// the p=11 generator orders (8,4).  The small representative below is
// therefore enough to identify their common reduction without constructing
// enormous exact coordinates for the five high-height rows.
Analyze("V2",[Q!-9,14,49],Q!-2023/162,Q!3/17,<1,3>,
        [<3,-1,1>],[11,13]);
UnsetLogFile();
quit;
