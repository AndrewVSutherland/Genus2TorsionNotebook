//////////////////////////////////////////////////////////////////////
// Cheap large-prime filter for the six V2 MW classes surviving p<=97.
// This avoids constructing enormous exact rational coordinates.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned FamilyName then FamilyName:="V2"; end if;
if not assigned input_file then
    input_file:=Sprintf("results/target_22224_%o_rank2_mwsieve_13_N1001_modular.tsv",FamilyName);
end if;
if not assigned log_file then
    log_file:=Sprintf("results/target_22224_%o_survivor_prime_filter.log",FamilyName);
end if;
if not assigned PrimeBound then PrimeBound:=997;
elif Type(PrimeBound) eq MonStgElt then PrimeBound:=StringToInteger(PrimeBound); end if;
SetLogFile(log_file:Overwrite:=true);

if FamilyName eq "V2" then
    fixed:=[Q!-9,Q!14,Q!49]; d0:=Q!-2023/162;
    qi:=1; qj:=3; remaining:=2;
elif FamilyName eq "W2" then
    fixed:=[Q!-8,Q!9,Q!25]; d0:=Q!-2209/338;
    qi:=1; qj:=2; remaining:=3;
elif FamilyName eq "AA1" then
    fixed:=[Q!-64,Q!-40,Q!65]; d0:=Q!9800/117;
    qi:=1; qj:=3; remaining:=2;
else
    error "unsupported family";
end if;
Rx:=(fixed[qi]+d0*T^2)/(fixed[qi]+d0);
Ry:=(fixed[qj]+d0*T^2)/(fixed[qj]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
// Match the bases and torsion-coset order printed by the originating sieves.
if FamilyName eq "V2" then
    W:=E![Q!-67,Q!0,Q!1];
    G1:=E![Q!573/4,Q!17661/8,Q!1];
    G2:=E![Q!-3,Q!888,Q!1];
    tors:=[E!0,W];
elif FamilyName eq "W2" then
    W:=E![Q!33,Q!272,Q!1]; assert Order(W) eq 4;
    G1:=E![Q!288,Q!4913,Q!1];
    G2:=E![Q!577,Q!13872,Q!1];
    tors:=[E!0,W,2*W,3*W];
else
    W:=E![Q!-45/4,Q!41/8,Q!1]; assert Order(W) eq 2;
    G1:=E![Q!1045,Q!75527,Q!1];
    G2:=E![Q!5,Q!-8453,Q!1];
    tors:=[E!0,W];
end if;
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);

function ReduceRat(q,F)
    q:=Q!q; ell:=Z!Characteristic(F);
    if Denominator(q) mod ell eq 0 then return false,F!0; end if;
    return true,F!Numerator(q)/F!Denominator(q);
end function;
function ReducePoint(PQ,EF)
    if PQ eq Parent(PQ)!0 then return true,EF!0; end if;
    cs:=[Q!PQ[i]:i in [1..3]]; den:=LCM([Denominator(z):z in cs]);
    ints:=[Z!(z*den):z in cs]; g:=GCD([Abs(z):z in ints]);
    if g ne 0 then ints:=[z div g:z in ints]; end if;
    try return true,EF![BaseRing(EF)!z:z in ints];
    catch e return false,EF!0; end try;
end function;
function ReduceMap(polys,F)
    out:=[];
    try for h in polys do Append(~out,ChangeRing(h,F)); end for;
    catch e return false,[]; end try;
    return true,out;
end function;
function EvalMap(polys,coords)
    vals:=[Evaluate(h,coords):h in polys];
    if &and[x eq 0:x in vals] then return false,vals; end if;
    return true,vals;
end function;
function TFromPoint(PF,rawF,curveF,F)
    if PF eq Parent(PF)!0 then return true,false,F!1; end if;
    ok,raw:=EvalMap(rawF,[PF[i]:i in [1..3]]);
    if not ok then return true,true,F!0; end if;
    ok,cp:=EvalMap(curveF,raw);
    if not ok or cp[3] eq 0 then return true,true,F!0; end if;
    return true,false,cp[1]/cp[3];
end function;

rows:=[];
for line in Split(Read(input_file),"\n") do
    if #line eq 0 or line[1] notin {"-","0","1","2","3","4","5","6","7","8","9"} then continue; end if;
    z:=Split(line,"\t"); if #z lt 3 then continue; end if;
    Append(~rows,<StringToInteger(z[1]),StringToInteger(z[2]),StringToInteger(z[3])>);
end for;
print "SURVIVOR_PRIME_FILTER_START","family",FamilyName,"rows",rows,"bound",PrimeBound;

alive:=[true:i in [1..#rows]]; killed:=0;
for ell in PrimesInInterval(101,PrimeBound) do
    if &and[not z:z in alive] then break; end if;
    F:=GF(ell); ainv:=[]; good:=true;
    for q in aInvariants(E) do
        ok,v:=ReduceRat(q,F); if not ok then good:=false; break; end if;
        Append(~ainv,v);
    end for;
    if not good then continue; end if;
    try EF:=EllipticCurve(ainv); catch e continue; end try;
    ok1,g1:=ReducePoint(G1,EF); ok2,g2:=ReducePoint(G2,EF);
    if not ok1 or not ok2 then continue; end if;
    torsF:=[];
    for w in tors do ok,wp:=ReducePoint(w,EF); if not ok then good:=false; break; end if; Append(~torsF,wp); end for;
    if not good then continue; end if;
    okr,rawF:=ReduceMap(rawMap,F); okc,curveF:=ReduceMap(curveMap,F);
    if not okr or not okc then continue; end if;
    P<X>:=PolynomialRing(F);
    for i in [1..#rows] do if alive[i] then
        m:=rows[i][1]; n:=rows[i][2]; ti:=rows[i][3];
        ep:=m*g1+n*g2+torsF[ti]; ok,isinf,tt:=TFromPoint(ep,rawF,curveF,F);
        if not ok or isinf then continue; end if;
        den:=F!(fixed[remaining]+d0);
        if den eq 0 then continue; end if;
        r3:=(F!fixed[remaining]+F!d0*tt^2)/den;
        bad:=not IsSquare(r3);
        if not bad then
            vals:=[F!z:z in fixed] cat [F!d0*tt^2];
            f:=X*&*[X+z^2:z in vals];
            if Degree(f) eq 5 and Discriminant(f) ne 0 and r3 ne 0 then
                nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
                bad:=nj mod 3 ne 0;
            end if;
        end if;
        if bad then
            alive[i]:=false; killed+:=1;
            print "CLASS_KILLED","family",FamilyName,"row",rows[i],"prime",ell,"t",tt;
        end if;
    end if; end for;
end for;
surv:=[rows[i]:i in [1..#rows]|alive[i]];
print "SURVIVOR_PRIME_FILTER_DONE","family",FamilyName,"killed",killed,"survivors",surv;
for row in surv do
    ep:=row[1]*G1+row[2]*G2+tors[row[3]];
    try cp:=Einv(minmapinv(ep));
        if cp[3] ne 0 then
            tt:=Q!(cp[1]/cp[3]); dd:=d0*tt^2;
            sq:=[IsSquare((z+dd)/(z+d0)):z in fixed];
            print "EXACT_SURVIVOR","family",FamilyName,"row",row,"t",tt,"square_flags",sq;
        end if;
    catch e
        print "EXACT_SURVIVOR_ERROR","family",FamilyName,"row",row,e;
    end try;
end for;
UnsetLogFile(); quit;
