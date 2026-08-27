//////////////////////////////////////////////////////////////////////
// Auxiliary-prime filter for the V1 N=10000 rows that genuinely meet a
// p=13 contact+full-cover component.  A nonsquare third quotient or a good
// reduction with 3 not dividing #J(F_l) is a rigorous global obstruction.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(8);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
input_file:="results/target_22224_V1_MW_p13_components.tsv";
output_file:="results/target_22224_V1_component_auxprime_filter.tsv";
log_file:="results/target_22224_V1_component_auxprime_filter.log";
SetLogFile(log_file:Overwrite:=true);

fixed:=[Q!-18,Q!20,Q!75]; d0:=Q!-1470/121;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
GQ:=[E![Q!3603,Q!216600,Q!1],E![Q!93,Q!-2100,Q!1]];
torsQ:=[E!0,E![Q!-7,Q!0,Q!1]];
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
    catch err return false,EF!0; end try;
end function;
function ReduceMap(polys,F)
    out:=[];
    try
        for h in polys do Append(~out,ChangeRing(h,F)); end for;
        return true,out;
    catch err
        allcoeffs:=&cat[Coefficients(h):h in polys];
        den:=LCM([Denominator(Q!c):c in allcoeffs]);
        scaled:=[den*h:h in polys];
        coeffs:=[Z!c:c in &cat[Coefficients(h):h in scaled]];
        cont:=GCD([Abs(c):c in coeffs]);
        if cont eq 0 then return false,[]; end if;
        scaled:=[h/cont:h in scaled]; out:=[];
        try
            for h in scaled do Append(~out,ChangeRing(h,F)); end for;
        catch err2
            return false,[];
        end try;
        return true,out;
    end try;
end function;
function EvalMap(polys,coords)
    vals:=[Evaluate(h,coords):h in polys];
    if &and[z eq 0:z in vals] then return false,vals; end if;
    return true,vals;
end function;
function TFromPoint(PF,rawF,curveF,F)
    if PF eq Parent(PF)!0 then return true,false,F!1; end if;
    ok,raw:=EvalMap(rawF,[PF[i]:i in [1..3]]);
    if not ok then return false,true,F!0; end if;
    ok,cp:=EvalMap(curveF,raw);
    if not ok or cp[3] eq 0 then return false,true,F!0; end if;
    return true,false,cp[1]/cp[3];
end function;

// Read only genuine p=13 component hits.
rows:=[]; lines:=Split(Read(input_file),"\n");
for i in [2..#lines] do
    if #lines[i] eq 0 then continue; end if;
    w:=Split(lines[i],"\t"); if #w lt 9 or w[9] eq "none" then continue; end if;
    Append(~rows,<StringToInteger(w[1]),StringToInteger(w[2]),
                  StringToInteger(w[3]),w[9]>);
end for;
print "V1_AUX_FILTER_START","rows",rows;

out:=Open(output_file,"w");
fprintf out,"m\tn\ttorsion_coset\tcomponent\tstatus\tprime\treason\tT_mod_prime\tJ_order\n";
survivors:=[];
for row in rows do
    m:=row[1]; n:=row[2]; ti:=row[3]; component:=row[4]; killed:=false;
    tested:=0;
    if Max(Abs(m),Abs(n)) le 20 then
        try
            cpQ:=Einv(minmapinv(m*GQ[1]+n*GQ[2]+torsQ[ti]));
            if cpQ[3] ne 0 then
                ttQ:=Q!(cpQ[1]/cpQ[3]);
                r3Q:=(fixed[3]+d0*ttQ^2)/(fixed[3]+d0);
                if not IsSquare(r3Q) then
                    fprintf out,"%o\t%o\t%o\t%o\tkilled\t0\texact_third_nonsquare\t-1\t-1\n",
                            m,n,ti,component;
                    print "V1_AUX_KILL",row,"exact_third_nonsquare","T",ttQ;
                    killed:=true;
                end if;
            end if;
        catch err
            killed:=false;
        end try;
    end if;
    if killed then continue; end if;
    for ell in PrimesInInterval(101,2003) do
        if ell in {11,13} then continue; end if;
        F:=GF(ell); ainv:=[]; good:=true;
        for q in aInvariants(E) do ok,z:=ReduceRat(q,F);
            if not ok then good:=false; break; end if; Append(~ainv,z);
        end for;
        if not good then continue; end if;
        try EF:=EllipticCurve(ainv); catch err continue; end try;
        ok1,g1:=ReducePoint(GQ[1],EF); ok2,g2:=ReducePoint(GQ[2],EF);
        okt,tq:=ReducePoint(torsQ[ti],EF);
        okr,rawF:=ReduceMap(rawMap,F); okc,curveF:=ReduceMap(curveMap,F);
        if not (ok1 and ok2 and okt and okr and okc) then continue; end if;
        ep:=m*g1+n*g2+tq; okT,isinf,tt:=TFromPoint(ep,rawF,curveF,F);
        if not okT or isinf then continue; end if;
        tested+:=1;
        den:=F!(fixed[3]+d0);
        if den eq 0 then continue; end if;
        r3:=(F!fixed[3]+F!d0*tt^2)/den;
        if not IsSquare(r3) then
            fprintf out,"%o\t%o\t%o\t%o\tkilled\t%o\tthird_nonsquare\t%o\t-1\n",
                    m,n,ti,component,ell,Z!tt;
            print "V1_AUX_KILL",row,"p",ell,"third_nonsquare","T",Z!tt;
            killed:=true; break;
        end if;
        P<x>:=PolynomialRing(F);
        vals:=[F!fixed[1],F!fixed[2],F!fixed[3],F!d0*tt^2];
        f:=x*&*[x+z^2:z in vals];
        if Degree(f) ne 5 or Discriminant(f) eq 0 then continue; end if;
        nJ:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
        if nJ mod 3 ne 0 then
            fprintf out,"%o\t%o\t%o\t%o\tkilled\t%o\tJ_order_not_div3\t%o\t%o\n",
                    m,n,ti,component,ell,Z!tt,nJ;
            print "V1_AUX_KILL",row,"p",ell,"Jorder",nJ,"T",Z!tt;
            killed:=true; break;
        end if;
    end for;
    if not killed then
        Append(~survivors,row);
        fprintf out,"%o\t%o\t%o\t%o\tsurvives\t-1\tprimes_through_2003\t-1\t-1\n",
                m,n,ti,component;
        print "V1_AUX_SURVIVES",row,"tested",tested;
    end if;
end for;
delete out;
print "V1_AUX_FILTER_DONE","input_rows",#rows,"survivors",survivors,
      "output",output_file;
UnsetLogFile(); quit;
