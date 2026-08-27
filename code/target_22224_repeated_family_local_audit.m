//////////////////////////////////////////////////////////////////////
// Component-complete first-pass local audit for every repeated
// A(2,2,2,8) fiber found in the complete B=10000 transverse box.
//
// At every coefficient-good prime p != 2,3 we enumerate all integral
// T residues satisfying the four full-cover square conditions.  A smooth
// residue is killed when 3 does not divide #J(F_p).  Singular residues and
// residues with 3 | #J(F_p) are retained.  The v_p(T)<0 chart is killed
// only when its four leading cover units are not all squares.
//
// Thus DEAD means a rigorous local obstruction for the whole rational
// one-parameter family; LIVE merely means this inexpensive audit does not
// decide the family.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Q:=Rationals(); Z:=Integers();

if not assigned PrimeList then
    PrimeList:=[5,7,11,13,17,19,23,29,31,37,41,43];
elif Type(PrimeList) eq MonStgElt then
    PrimeList:=[StringToInteger(s):s in Split(PrimeList,",")];
end if;
if not assigned FamilyNames then
    FamilyNames:="";
end if;
selectedNames:=Type(FamilyNames) eq MonStgElt and #FamilyNames gt 0
               select SequenceToSet(Split(FamilyNames,",")) else {};
if not assigned log_file then
    log_file:="results/target_22224_repeated_family_local_audit.log";
end if;
if not assigned output_file then
    output_file:="results/target_22224_repeated_family_local_audit.tsv";
end if;
SetLogFile(log_file:Overwrite:=true);

fibers := [
 <"F1",[Q!-1470,-630,336],Q!25,Q!5>,
 <"F2",[Q!-720,20,300],Q!-363,Q!17/11>,
 <"F3",[Q!-612,34,289],Q!-338,Q!25/13>,
 <"F4",[Q!-126,28,49],Q!-50,Q!17/5>,
 <"F5",[Q!-112,14,49],Q!-50,Q!13/5>,
 <"F6",[Q!-50,30,45],Q!-48,Q!2>,
 <"F7",[Q!-18,1,16],Q!-50,Q!29/5>,
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
 <"N11",[Q!-64,120,169],Q!-1352/15,Q!1/13>,
 <"N12",[Q!-28,50,55],Q!-3610/77,Q!1/19>,
 <"O1",[Q!-1,4,6],Q!-98/27,Q!3/7>,
 <"O2",[Q!-1,6,9],Q!-507/98,Q!7/13>,
 <"O3",[Q!-3,-2,12],Q!162/49,Q!7/3>,
 <"O4",[Q!-9,10,15],Q!-75/8,Q!1/2>,
 <"O5",[Q!-15,-10,16],Q!50/3,Q!5/3>,
 <"O6",[Q!-15,-1,25],Q!605/27,Q!21/11>,
 <"O7",[Q!-9,-5,30],Q!75/8,Q!5/2>,
 <"O8",[Q!-22,55,70],Q!-1372/25,Q!5/7>,
 <"O9",[Q!1,55,99],Q!5/9,Q!15>,
 <"M1",[Q!-16,-5,30],Q!50/3,Q!11/3>,
 <"M2",[Q!-1,25,40],Q!-605/32,Q!17/22>,
 <"M3",[Q!-25,-1,65],Q!405/13,Q!37/27>,
 <"M4",[Q!-49,-2,100],Q!3362/49,Q!61/41>,
 <"M5",[Q!99,244,4026],Q!6,Q!29>,
 <"R1",[Q!-1,2,25],Q!-529/338,Q!13/23>,
 <"R2",[Q!-44,-16,49],Q!676/11,Q!43/13>,
 <"R3",[Q!-17,-16,50],Q!338/17,Q!53/13>,
 <"R4",[Q!-7,16,56],Q!-392/25,Q!5/13>,
 <"R5",[Q!-23,-18,64],Q!578/23,Q!65/17>,
 <"R6",[Q!-3,5,75],Q!-405/121,Q!11/21>,
 <"R7",[Q!-47,-4,98],Q!2738/47,Q!65/37>,
 <"S1",[Q!-4,9,30],Q!-2166/245,Q!7/19>,
 <"S3",[Q!-9,35,65],Q!-3025/91,Q!5/11>,
 <"U1",[Q!-119,124,126],Q!-1054/9,Q!3/7>,
 <"U2",[Q!-1071,-1054,1116],Q!1134,Q!7/3>,
 <"V1",[Q!-18,20,75],Q!-1470/121,Q!11/49>,
 <"V2",[Q!-9,14,49],Q!-2023/162,Q!3/17>,
 <"W1",[Q!-2254,-2162,2303],Q!4900,Q!7/5>,
 <"W2",[Q!-8,9,25],Q!-2209/338,Q!13/47>,
 <"X1",[Q!-49,-21,81],Q!6069/121,Q!143/17>,
 <"Y1",[Q!-15,28,50],Q!-5290/189,Q!3/23>,
 <"Z1",[Q!-169,-120,225],Q!5415/32,Q!26/19>,
 <"AA1",[Q!-64,-40,65],Q!9800/117,Q!9/7>
];

function ReduceRat(q,F)
    q:=Q!q; p:=Z!Characteristic(F);
    if Denominator(q) mod p eq 0 then return false,F!0; end if;
    return true,F!Numerator(q)/F!Denominator(q);
end function;

function CoverRadicands(vals)
    a,b,c,d:=Explode(vals);
    return [a*b*c*d,
            a*(a+b)*(a+c)*(a+d),
            b*(b+a)*(b+c)*(b+d),
            c*(c+a)*(c+b)*(c+d)];
end function;

out:=Open(output_file,"w");
fprintf out,"family\tp\tstatus\tfinite_cover\tfinite_smooth\tfinite_live\tfinite_singular\tinfinity_full\tlive_T\tJ_orders\n";
deadNames:={}; liveNames:={}; skippedNames:={};

print "REPEATED_FAMILY_LOCAL_AUDIT_START","families",#fibers,"primes",PrimeList;
for fam in fibers do
    name:=fam[1]; fixed:=fam[2]; d0:=fam[3]; knownT:=fam[4];
    if #selectedNames gt 0 and name notin selectedNames then continue; end if;
    famDead:=false; usable:=0;
    for p in PrimeList do
        if p in {2,3} then continue; end if;
        F:=GF(p); P<X>:=PolynomialRing(F);
        ok:=true; ff:=[];
        for z in fixed do
            good,zz:=ReduceRat(z,F); ok and:=good and zz ne 0; Append(~ff,zz);
        end for;
        good,dd:=ReduceRat(d0,F); ok and:=good and dd ne 0;
        if not ok then
            fprintf out,"%o\t%o\tSKIP_COEFFICIENT\t-1\t-1\t-1\t-1\t-1\t[]\t[]\n",name,p;
            continue;
        end if;
        usable+:=1;
        finiteCover:=0; finiteSmooth:=0; finiteSing:=0;
        liveT:=[]; jorders:=[];
        for tt in F do
            vals:=ff cat [dd*tt^2];
            rr:=CoverRadicands(vals);
            if not &and[IsSquare(r):r in rr] then continue; end if;
            finiteCover+:=1;
            f:=X*&*[X+z^2:z in vals];
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                finiteSing+:=1; Append(~liveT,Z!tt); Append(~jorders,-1);
                continue;
            end if;
            finiteSmooth+:=1;
            nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
            if nj mod 3 eq 0 then Append(~liveT,Z!tt); Append(~jorders,nj); end if;
        end for;
        A,B,C:=Explode(ff); D:=dd;
        infUnits:=[A*B*C*D,
                   A*(A+B)*(A+C)*D,
                   B*(B+A)*(B+C)*D,
                   C*(C+A)*(C+B)*D];
        infFull:=&and[IsSquare(r):r in infUnits];
        if infFull then Append(~liveT,-1); Append(~jorders,-2); end if;
        status:=(#liveT eq 0) select "DEAD" else "LIVE";
        fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                name,p,status,finiteCover,finiteSmooth,#liveT,finiteSing,
                infFull select 1 else 0,liveT,jorders;
        if status eq "DEAD" then
            famDead:=true; Include(~deadNames,name);
            print "FAMILY_DEAD",name,"p",p,"finite_cover",finiteCover,
                  "finite_smooth",finiteSmooth,"infinity_full",infFull;
            break;
        end if;
    end for;
    if not famDead then
        if usable eq 0 then Include(~skippedNames,name); else Include(~liveNames,name); end if;
        print "FAMILY_NOT_KILLED",name,"usable_primes",usable,"known_T",knownT;
    end if;
end for;
delete out;
print "REPEATED_FAMILY_LOCAL_AUDIT_DONE","dead",#deadNames,Sort(Setseq(deadNames)),
      "not_killed",#liveNames,Sort(Setseq(liveNames)),
      "skipped",#skippedNames,Sort(Setseq(skippedNames)),"output",output_file;
UnsetLogFile(); quit;
