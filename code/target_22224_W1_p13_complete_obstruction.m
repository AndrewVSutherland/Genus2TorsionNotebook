//////////////////////////////////////////////////////////////////////
// Component-complete p=13 obstruction for the repeated fiber
//
//   W1 = (-2254,-2162,2303,4900*T^2).
//
// The four full-cover radicands are tested on every finite residue of T.
// For every surviving smooth residue we compute #J(F_13) and also enumerate
// the corrected affine cubic-contact equations with M=L^2, L nonzero.
// The nonintegral T chart is tested after z=1/T and projective rescaling.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Q:=Rationals(); Z:=Integers(); p:=13; F:=GF(p); P<X>:=PolynomialRing(F);
fixed:=[Z!-2254,Z!-2162,Z!2303]; d0:=Z!4900;
log_file:="results/target_22224_W1_p13_complete_obstruction.log";
output_file:="results/target_22224_W1_p13_complete_obstruction.tsv";
bank_file_pre:="results/target_22224_W1_rank3_periodic_N1000.tsv";
bank_file:="results/target_22224_W1_rank3_periodic_N1000_p167.tsv";
bank_output:="results/target_22224_W1_rank3_periodic_N1000_p13complete.tsv";
coefficient_mask:="results/target_22224_W1_p13_coefficient_mask.tsv";
SetLogFile(log_file:Overwrite:=true);

function CoverRadicands(vals)
    a,b,c,d:=Explode(vals);
    return [a*b*c*d,
            a*(a+b)*(a+c)*(a+d),
            b*(b+a)*(b+c)*(b+d),
            c*(c+a)*(c+b)*(c+d)];
end function;

function Elementary(vals)
    a,b,c,d:=Explode([z^2:z in vals]);
    return a+b+c+d,
           a*b+a*c+a*d+b*c+b*d+c*d,
           a*b*c+a*b*d+a*c*d+b*c*d,
           a*b*c*d;
end function;

function ContactCount(vals)
    e1,e2,e3,e4:=Elementary(vals); count:=0; samples:=[];
    for L in F do if L ne 0 then
        M:=L^2;
        for U in F do for v in F do
            PP:=4*M*e1+12*(U^2+v^2)-(M+3*U)^2;
            G1:=(M+3*U)*PP+16*v^3-8*U^3-48*U*v^2-8*M*e2;
            G2:=PP^2+64*(M+3*U)*v^3-192*(U^2*v^2+v^4)-64*M*e3;
            G3:=PP*v^3-12*U*v^4-4*M*e4;
            if G1 eq 0 and G2 eq 0 and G3 eq 0 then
                count+:=1; Append(~samples,<Z!L,Z!U,Z!v>);
            end if;
        end for; end for;
    end if; end for;
    return count,samples;
end function;

print "W1_P13_COMPLETE_START","fiber",fixed cat [d0],"p",p;
out:=Open(output_file,"w");
fprintf out,"chart\tT_residue\tcover_radicands\tsmooth\tJ_order\tJ_order_mod3\taffine_contact_count\n";

finiteCover:=[]; smoothCover:=0; totalContact:=0;
for tt in F do
    vals:=[F!z:z in fixed] cat [F!d0*tt^2];
    rr:=CoverRadicands(vals);
    if not &and[IsSquare(r):r in rr] then continue; end if;
    Append(~finiteCover,Z!tt);
    f:=X*&*[X+z^2:z in vals];
    smooth:=Degree(f) eq 5 and Discriminant(f) ne 0;
    nj:=-1;
    if smooth then
        smoothCover+:=1;
        C:=HyperellipticCurve(f); nj:=Z!Evaluate(LPolynomial(C),1);
    end if;
    nc,samples:=ContactCount(vals); totalContact+:=nc;
    fprintf out,"finite\t%o\t%o\t%o\t%o\t%o\t%o\n",
            Z!tt,[Z!r:r in rr],smooth select 1 else 0,nj,
            nj eq -1 select -1 else nj mod 3,nc;
    print "FINITE_COMPONENT","T",Z!tt,"signed_vals",[Z!z:z in vals],
          "radicands",[Z!r:r in rr],"smooth",smooth,
          "J_order",nj,"mod3",nj eq -1 select -1 else nj mod 3,
          "contact_count",nc,"samples",samples;
end for;

// For v_13(T)<0 put z=1/T and rescale the signed branch tuple by z^2:
// [A*z^2,B*z^2,C*z^2,D].  Divide the four cover radicands by z^6.
A,B,C:=Explode([F!z:z in fixed]); D:=F!d0;
infUnits:=[A*B*C*D,
           A*(A+B)*(A+C)*D,
           B*(B+A)*(B+C)*D,
           C*(C+A)*(C+B)*D];
infFull:=&and[IsSquare(r):r in infUnits];
fprintf out,"infinity\t-1\t%o\t-1\t-1\t-1\t-1\n",[Z!r:r in infUnits];
print "INFINITY_COMPONENT","units",[Z!r:r in infUnits],
      "square_flags",[IsSquare(r):r in infUnits],"full_cover",infFull;

// The p=13 obstruction is intrinsic, so every coefficient class is killed.
// Audit the strongest existing periodic survivor file and emit an empty bank.
prelines:=Split(Read(bank_file_pre),"\n"); preRows:=0;
for i in [2..#prelines] do if #prelines[i] gt 0 then preRows+:=1; end if; end for;
lines:=Split(Read(bank_file),"\n"); bankRows:=0;
for i in [2..#lines] do if #lines[i] gt 0 then bankRows+:=1; end if; end for;
bout:=Open(bank_output,"w");
fprintf bout,"m\tn\tk\ttorsion_coset\n"; delete bout;
cm:=Open(coefficient_mask,"w");
fprintf cm,"modulus\tcondition\tsurviving_coefficient_classes\treason\n";
fprintf cm,"13\tFALSE\t0\tno finite or infinity full-cover/contact component\n";
delete cm;

delete out;
assert Sort(finiteCover) eq [1,4,9,12];
assert smoothCover eq 4;
assert totalContact eq 0;
assert not infFull;
print "MW_CONGRUENCE_MASK","modulus",13,"allowed_classes","EMPTY",
      "reason","no finite or infinity full-cover/contact component",
      "output",coefficient_mask;
print "PERIODIC_BANK_AUDIT","N",1000,"raw_classes",2*2001^3,
      "initial_periodic_rows",preRows,"through_p167_rows",bankRows,
      "p13_component_survivors",0,"output",bank_output;
print "W1_P13_COMPLETE_DONE","finite_cover_residues",Sort(finiteCover),
      "smooth_cover",smoothCover,"affine_contact_total",totalContact,
      "infinity_full_cover",infFull,"global_W1_survivors",0,
      "output",output_file;
UnsetLogFile(); quit;
