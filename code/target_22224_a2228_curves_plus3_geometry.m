//////////////////////////////////////////////////////////////////////
// Correct cubic-contact pullbacks to six explicit P^1 curves on
// the full A(2,2,2,8) cover.
//
// For f=x^5+e1*x^4+e2*x^3+e3*x^2+e4*x and
//
//   q=x^2+U*x+v^2,
//   h=(1/L)x^3+N*x^2+R*x+v^3/L,
//
// the identity h^2-f=(1/L^2)q^3 gives M=L^2 and
//
//   P  = 4*M*e1+12*(U^2+v^2)-(M+3*U)^2,
//   N  = (M+3*U)/(2*L),       R=P/(8*L),
//   G1 = (M+3*U)*P+16*v^3-8*M*e2-8*U^3-48*U*v^2,
//   G2 = P^2+64*(M+3*U)*v^3-64*M*e3
//          -192*(U^2*v^2+v^4),
//   G3 = P*v^3-4*M*e4-12*U*v^4.
//
// This corrects the older Aaux formula, which doubled the
// -(M+3U)^2 term.  Mode summary writes the exact pulled-back equations;
// mode validate checks the formulas against #J(F_p) and constructs an
// order-three divisor for every contact witness at small primes.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then mode := "summary"; end if;
if not assigned output_file then
    output_file := "results/target_22224_a2228_curves_plus3_equations.txt";
end if;
if not assigned PrimeList then PrimeList := [11,13,17,19,23];
elif Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s):s in Split(PrimeList,",")];
end if;

Q := Rationals(); Z := Integers();

function FamilyName(family)
    return ["Filip1","Filip2","Adam","Filip3",
            "Product_t1_2","Product_t2_3"][family];
end function;

function FamilyTuple(family,t)
    K := Parent(t);
    if family eq 1 then
        return [-(t^2+t+1)^2,-4*t*(t+1)^2,4*t*(t+1),4*t^2*(t+1)];
    elif family eq 2 then
        den := t^4+2*t^3-t^2-2*t+1;
        return [-(t^4-2*t^3-t^2+2*t+1)/den,-1/t,K!1,t];
    elif family eq 3 then
        return [-t*(t+1)/(t-1),K!1,-t^2,t*(t-1)/(t+1)];
    elif family eq 4 then
        den := (t^2-2*t-1)*(t^2+1);
        return [-t*(t+1)^2*(t-1)/den,-t^2,K!1,t];
    elif family in {5,6} then
        A := family eq 5 select 16 else 144;
        B := family eq 5 select 9 else 25;
        s := (B-A*t)/(A-B*t);
        return [K!1,t,s,t*s];
    end if;
    error "family must lie in [1..6]";
end function;

function HomogeneousTuple(family,n,d)
    if family eq 1 then
        return [-(n^2+n*d+d^2)^2,
                -4*n*(n+d)^2*d,4*n*(n+d)*d^2,4*n^2*(n+d)*d];
    elif family eq 2 then
        N := n^4-2*n^3*d-n^2*d^2+2*n*d^3+d^4;
        D := n^4+2*n^3*d-n^2*d^2-2*n*d^3+d^4;
        return [-n*d*N,-d^2*D,n*d*D,n^2*D];
    elif family eq 3 then
        return [-n*(n+d)^2*d,(n^2-d^2)*d^2,
                -n^2*(n^2-d^2),n*(n-d)^2*d];
    elif family eq 4 then
        N := n*(n+d)^2*(n-d);
        D := (n^2-2*n*d-d^2)*(n^2+d^2);
        return [-N*d^2,-n^2*D,d^2*D,n*d*D];
    elif family in {5,6} then
        A := family eq 5 select 16 else 144;
        B := family eq 5 select 9 else 25;
        D := A*d-B*n; N := B*d-A*n;
        return [d*D,n*D,d*N,n*N];
    end if;
    error "family must lie in [1..6]";
end function;

function PrimitiveNumerator(g,R)
    h := R!Numerator(g);
    raw := [c:c in Coefficients(h)|c ne 0];
    if #raw eq 0 then return h; end if;
    den := LCM([Denominator(c):c in raw]);
    h *:= den;
    cc := [Z!c:c in Coefficients(h)|c ne 0];
    if #cc eq 0 then return h; end if;
    content := GCD([Abs(c):c in cc]);
    if content gt 1 then h /:= content; end if;
    return h;
end function;

function PulledBackEquations(family)
    R<t,L,U,v> := PolynomialRing(Q,4);
    K := FieldOfFractions(R); PX<x> := PolynomialRing(K);
    vals := FamilyTuple(family,K!t);
    f := x*&*[x+a^2:a in vals];
    e1,e2,e3,e4 := Explode([Coefficient(f,i):i in [4,3,2,1]]);
    M := L^2;
    PP := 4*M*e1+12*(U^2+v^2)-(M+3*U)^2;
    G1 := (M+3*U)*PP+16*v^3-8*M*e2-8*U^3-48*U*v^2;
    G2 := PP^2+64*(M+3*U)*v^3-64*M*e3
          -192*(U^2*v^2+v^4);
    G3 := PP*v^3-4*M*e4-12*U*v^4;
    return [PrimitiveNumerator(g,R):g in [G1,G2,G3]],f;
end function;

function BadReductionFactorization(family)
    K<t> := RationalFunctionField(Q);
    vals := FamilyTuple(family,t);
    bad := &*vals;
    for i in [1..4] do for j in [1..i-1] do
        bad *:= (vals[i]-vals[j])*(vals[i]+vals[j]);
    end for; end for;
    return Factorization(Numerator(bad)),Factorization(Denominator(bad));
end function;

function SmoothTuple(vals)
    sq := [z^2:z in vals];
    return &and[z ne 0:z in sq] and #Set(sq) eq 4;
end function;

function ContactWitness(vals)
    F := Universe(vals); PX<x> := PolynomialRing(F);
    f := x*&*[x+a^2:a in vals];
    e1,e2,e3,e4 := Explode([Coefficient(f,i):i in [4,3,2,1]]);
    for L in F do
      if L eq 0 then continue; end if;
      M := L^2;
      for U in F do for v in F do
        PP := 4*M*e1+12*(U^2+v^2)-(M+3*U)^2;
        G1 := (M+3*U)*PP+16*v^3-8*M*e2-8*U^3-48*U*v^2;
        if G1 ne 0 then continue; end if;
        G2 := PP^2+64*(M+3*U)*v^3-64*M*e3
              -192*(U^2*v^2+v^4);
        if G2 ne 0 then continue; end if;
        G3 := PP*v^3-4*M*e4-12*U*v^4;
        if G3 ne 0 then continue; end if;
        q := x^2+U*x+v^2;
        if Degree(GCD(q,f)) ne 0 then continue; end if;
        N := (M+3*U)/(2*L); Rr := PP/(8*L); S := v^3/L;
        h := (1/L)*x^3+N*x^2+Rr*x+S;
        assert h^2-f eq (1/M)*q^3;
        C := HyperellipticCurve(f); J := Jacobian(C); D := J![q,h mod q];
        assert 3*D eq J!0 and D ne J!0;
        return true,<L,U,v,q,h>;
      end for; end for;
    end for;
    return false,<F!0,F!0,F!0,PX!0,PX!0>;
end function;

if mode eq "summary" then
    out := Open(output_file,"w");
    for family in [1..6] do
        Fs,f := PulledBackEquations(family);
        badnum,badden := BadReductionFactorization(family);
        fprintf out,"FAMILY %o %o\n",family,FamilyName(family);
        fprintf out,"curve = %o\n",f;
        fprintf out,"bad_numerator_factorization = %o\n",badnum;
        fprintf out,"bad_denominator_factorization = %o\n",badden;
        for i in [1..3] do
            g := Fs[i];
            fprintf out,"G%o = %o\n",i,g;
            fprintf out,"G%o_summary total_degree=%o degree_t=%o degree_L=%o degree_U=%o degree_v=%o terms=%o\n",
                i,TotalDegree(g),Degree(g,1),Degree(g,2),Degree(g,3),Degree(g,4),#Terms(g);
        end for;
        fprintf out,"\n";
        print "PULLBACK_SUMMARY",FamilyName(family),
              [<TotalDegree(g),Degree(g,1),Degree(g,2),Degree(g,3),Degree(g,4),#Terms(g)>:g in Fs];
        print "BAD_FACTORS",FamilyName(family),badnum,badden;
    end for;
    delete out;
    print "PULLBACK_EQUATIONS_FILE",output_file;
elif mode eq "validate" then
    for p in PrimeList do
        F := GF(p); checked := 0; group_positive := 0; contact_positive := 0;
        for family in [1..6] do
            fam_checked := 0; fam_group := 0; fam_contact := 0;
            for t in F do
                vals := HomogeneousTuple(family,t,F!1);
                if not SmoothTuple(vals) then continue; end if;
                fam_checked +:= 1; checked +:= 1;
                PX<x> := PolynomialRing(F);
                f := x*&*[x+a^2:a in vals];
                has_group := (Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1)) mod 3 eq 0;
                has_contact,w := ContactWitness(vals);
                assert has_group eq has_contact;
                if has_group then fam_group +:= 1; group_positive +:= 1; end if;
                if has_contact then fam_contact +:= 1; contact_positive +:= 1; end if;
            end for;
            print "CONTACT_VALIDATE","p",p,"family",family,
                  "checked",fam_checked,"group_positive",fam_group,
                  "contact_positive",fam_contact;
        end for;
        print "CONTACT_VALIDATE_PRIME_DONE",p,checked,group_positive,contact_positive;
    end for;
else
    error "mode must be summary or validate";
end if;

quit;
