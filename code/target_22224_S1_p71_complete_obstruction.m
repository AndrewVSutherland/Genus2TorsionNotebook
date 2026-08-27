//////////////////////////////////////////////////////////////////////
// Component-complete p=71 obstruction for
//
//   S1 = (-4,9,30,(-2166/245)*T^2).
//
// We use the integral projectively equivalent signed tuple
//
//   (-980,2205,7350,-2166*T^2).
//
// Every finite full-cover residue is smooth.  Its genus-two Jacobian order
// is computed directly.  The nonintegral T chart is checked after z=1/T.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Z:=Integers(); p:=71; F:=GF(p); P<X>:=PolynomialRing(F);
fixed:=[Z!-980,Z!2205,Z!7350]; d0:=Z!-2166;
log_file:="results/target_22224_S1_p71_complete_obstruction.log";
output_file:="results/target_22224_S1_p71_complete_obstruction.tsv";
coefficient_mask:="results/target_22224_S1_p71_coefficient_mask.tsv";
bank_file:="results/target_22224_S1_rank1_modular_N100000000.tsv";
bank_output:="results/target_22224_S1_rank1_N100000000_p71complete.tsv";
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

// Independent audit of the corrected affine contact equations.  Since the
// curve is smooth and #J(F_71) is prime to 3, every genuine count must be 0.
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
                count+:=1;
                if #samples lt 20 then Append(~samples,<Z!L,Z!U,Z!v>); end if;
            end if;
        end for; end for;
    end if; end for;
    return count,samples;
end function;

print "S1_P71_COMPLETE_START","scaled_fiber",fixed cat [d0],"p",p;
out:=Open(output_file,"w");
fprintf out,"chart\tT_residue\tcover_radicands\tsigned_branches\tbranch_squares\tsmooth\tJ_order\tJ_order_mod3\taffine_contact_count\n";

finiteCover:=[]; smoothCover:=0; j3Cover:=[]; totalContact:=0;
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
        nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
        if nj mod 3 eq 0 then Append(~j3Cover,<Z!tt,nj>); end if;
    end if;
    nc,samples:=ContactCount(vals); totalContact+:=nc;
    fprintf out,"finite\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            Z!tt,[Z!r:r in rr],[Z!z:z in vals],[Z!(z^2):z in vals],
            smooth select 1 else 0,nj,
            nj eq -1 select -1 else nj mod 3,nc;
    print "FINITE_COMPONENT","T",Z!tt,"signed_vals",[Z!z:z in vals],
          "branch_squares",[Z!(z^2):z in vals],
          "radicands",[Z!r:r in rr],"smooth",smooth,
          "J_order",nj,"mod3",nj eq -1 select -1 else nj mod 3,
          "contact_count",nc,"contact_samples",samples;
end for;

// For v_71(T)<0 put z=1/T and rescale to [A*z^2,B*z^2,C*z^2,D].
// Removing the common z^6 gives these four limiting cover units.
A,B,C:=Explode([F!z:z in fixed]); D:=F!d0;
infUnits:=[A*B*C*D,
           A*(A+B)*(A+C)*D,
           B*(B+A)*(B+C)*D,
           C*(C+A)*(C+B)*D];
infFlags:=[IsSquare(r):r in infUnits]; infFull:=&and infFlags;
fprintf out,"infinity\t-1\t%o\t[]\t[]\t-1\t-1\t-1\t-1\n",
        [Z!r:r in infUnits];
print "INFINITY_COMPONENT","units",[Z!r:r in infUnits],
      "square_flags",infFlags,"full_cover",infFull;

// The obstruction lies on the target genus-two curve, below the quotient;
// hence it forbids every Mordell--Weil coefficient and torsion coset.
lines:=Split(Read(bank_file),"\n"); bankRows:=0;
for i in [2..#lines] do if #lines[i] gt 0 then bankRows+:=1; end if; end for;
bout:=Open(bank_output,"w"); fprintf bout,"m\tn\ttorsion_coset\n"; delete bout;
cm:=Open(coefficient_mask,"w");
fprintf cm,"modulus\tcondition\tsurviving_coefficient_classes\treason\n";
fprintf cm,"71\tFALSE\t0\tno finite or infinity full-cover/contact component\n";
delete cm;
delete out;

expected:=[1,4,9,20,33,34,37,38,51,62,67,70];
assert Sort(finiteCover) eq expected;
assert smoothCover eq #expected;
assert #j3Cover eq 0;
assert totalContact eq 0;
assert not infFull;
print "MW_CONGRUENCE_MASK","modulus",71,"allowed_classes","EMPTY",
      "reason","no finite or infinity full-cover/contact component",
      "output",coefficient_mask;
print "PERIODIC_BANK_AUDIT","coefficient_bound",100000000,
      "raw_classes",400000002,"input_rows",bankRows,
      "p71_component_survivors",0,"output",bank_output;
print "S1_P71_COMPLETE_DONE","finite_cover_residues",Sort(finiteCover),
      "smooth_cover",smoothCover,"J3_cover",j3Cover,
      "affine_contact_total",totalContact,
      "infinity_full_cover",infFull,"global_S1_survivors",0,
      "output",output_file;
UnsetLogFile(); quit;
