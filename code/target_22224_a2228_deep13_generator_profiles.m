//////////////////////////////////////////////////////////////////////
// Correct point-count profiles for the A(2,2,2,8)+3 search.
//
// The first version of the full-family cubic-contact mask had a factor
// two error in the x^4 coefficient comparison.  This file deliberately
// avoids that elimination: for each of four explicit rational curves on
// the full four-radicand A(2,2,2,8) cover, and each point of P^1(F_p), it
// records either bad branch reduction or the exact value #J(F_p) mod 3.
// A rational [2,2,2,24] specialization must be in an allowed row at every
// prime (bad reduction, or #J(F_p) divisible by 3).
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned PrimeList then
    PrimeList := [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71];
elif Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s):s in Split(PrimeList,",")];
end if;
if not assigned output_file then
    output_file := "results/target_22224_a2228_deep13_generator_profiles.tsv";
end if;
if not assigned log_file then
    log_file := "results/target_22224_a2228_deep13_generator_profiles.log";
end if;

SetLogFile(log_file : Overwrite := true);
Z := Integers();

// Homogeneous polynomial representatives.  With t=m/n these are common
// multiples of Filip 1, Filip 2, and the recovered Adam curve.  The formula
// called `filip3` in the older order-16 search is not generically on the full
// four-radicand cover, so it is intentionally excluded here.
function FamilyTuple(fam,m,n)
    if fam eq 1 then
        return [
            -(m^2+m*n+n^2)^2,
            -4*m*n*(m+n)^2,
            4*m*n^2*(m+n),
            4*m^2*n*(m+n)
        ];
    elif fam eq 2 then
        H := m^4+2*m^3*n-m^2*n^2-2*m*n^3+n^4;
        N := m^4-2*m^3*n-m^2*n^2+2*m*n^3+n^4;
        return [-N*m*n,-H*n^2,H*m*n,H*m^2];
    end if;
    // Family 3 is Adam:
    // [-t(t+1)/(t-1),1,-t^2,t(t-1)/(t+1)].
    return [
        -m*n*(m+n)^2,
        n^2*(m^2-n^2),
        -m^2*(m^2-n^2),
        m*n*(m-n)^2
    ];
end function;

function CoverRadicands(vals)
    a,b,c,d := Explode(vals);
    return [
        a*b*c*d,
        a*(a+b)*(a+c)*(a+d),
        b*(b+a)*(b+c)*(b+d),
        c*(c+a)*(c+b)*(c+d)
    ];
end function;

out := Open(output_file,"w");
fprintf out,"family\tp\tparameter\tboundary\tjorder\tjmod3\tallowed\ta\tb\tc\td\n";

for p in PrimeList do
    F := GF(p); P<x> := PolynomialRing(F);
    total := 0; boundary_count := 0; good_count := 0; allowed_count := 0;
    good_allowed := 0;
    for fam in [1..3] do
        fam_allowed := 0; fam_good_allowed := 0; fam_boundary := 0;
        params := [<F!r,F!1,IntegerToString(r)>:r in [0..p-1]];
        Append(~params,<F!1,F!0,"inf">);
        for rec in params do
            m,n,label := Explode(rec);
            vals := FamilyTuple(fam,m,n);
            total +:= 1;
            allzero := &and[z eq 0:z in vals];
            sq := [z^2:z in vals];
            boundary := allzero or (&or[z eq 0:z in vals]) or #Set(sq) ne 4;
            nord := 0; jmod := -1;
            if boundary then
                boundary_count +:= 1; fam_boundary +:= 1;
            else
                assert &and[IsSquare(z):z in CoverRadicands(vals)];
                f := x*&*[x+z^2:z in vals];
                assert Degree(f) eq 5 and Discriminant(f) ne 0;
                nord := Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
                jmod := nord mod 3;
                good_count +:= 1;
            end if;
            allowed := boundary or jmod eq 0;
            if allowed then
                allowed_count +:= 1; fam_allowed +:= 1;
                if not boundary then
                    good_allowed +:= 1; fam_good_allowed +:= 1;
                end if;
            end if;
            fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                fam,p,label,boundary select 1 else 0,nord,jmod,
                allowed select 1 else 0,Z!vals[1],Z!vals[2],Z!vals[3],Z!vals[4];
        end for;
        print "PROFILE_FAMILY","p",p,"family",fam,"boundary",fam_boundary,
              "allowed",fam_allowed,"good_allowed",fam_good_allowed;
    end for;
    print "PROFILE_PRIME","p",p,"total",total,"boundary",boundary_count,
          "good",good_count,"allowed",allowed_count,"good_allowed",good_allowed;
    if p in {11,13} then
        assert good_allowed eq 0;
    end if;
end for;

delete out;
print "TARGET_22224_A2228_PROFILE_DONE","file",output_file;
UnsetLogFile();
quit;
