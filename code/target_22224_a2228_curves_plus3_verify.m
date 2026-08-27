//////////////////////////////////////////////////////////////////////
// Exact verifier for the denominator-aware P^1 CRT sieve.
//
// It first extends the reduction test beyond the primes used to create the
// masks.  Every survivor is sent to TorsionSubgroup; an order-three Mumford
// class is then converted back to an exact corrected cubic-contact identity.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned input_file then
    input_file := "results/target_22224_a2228_curves_plus3_crt_h100000.tsv";
end if;
if not assigned log_file then
    log_file := "results/target_22224_a2228_curves_plus3_verify.log";
end if;
if not assigned max_exact then max_exact := 100;
elif Type(max_exact) eq MonStgElt then max_exact := StringToInteger(max_exact);
end if;

SetLogFile(log_file : Overwrite := true);
Q := Rationals(); Z := Integers(); P<x> := PolynomialRing(Q);

function FamilyTuple(family,t)
    if family eq 1 then
        return [-(t^2+t+1)^2,-4*t*(t+1)^2,4*t*(t+1),4*t^2*(t+1)];
    elif family eq 2 then
        den := t^4+2*t^3-t^2-2*t+1;
        return [-(t^4-2*t^3-t^2+2*t+1)/den,-1/t,Q!1,t];
    elif family eq 3 then
        return [-t*(t+1)/(t-1),Q!1,-t^2,t*(t-1)/(t+1)];
    elif family eq 4 then
        den := (t^2-2*t-1)*(t^2+1);
        return [-t*(t+1)^2*(t-1)/den,-t^2,Q!1,t];
    elif family in {5,6} then
        A := family eq 5 select 16 else 144;
        B := family eq 5 select 9 else 25;
        s := (B-A*t)/(A-B*t);
        return [Q!1,t,s,t*s];
    end if;
    error "bad family";
end function;

function Usable(vals)
    return #vals eq 4 and &and[z ne 0:z in vals]
           and #Set([z^2:z in vals]) eq 4;
end function;

function CurvePolynomial(vals)
    return x*&*[x+z^2:z in vals];
end function;

function ReadRows(filename)
    rows := []; lines := Split(Read(filename),"\n");
    for line in lines do
        if #line eq 0 or line[1] notin "123456" then continue; end if;
        c := Split(line,"\t");
        if #c ge 3 then
            Append(~rows,<StringToInteger(c[1]),StringToInteger(c[2]),
                          StringToInteger(c[3])>);
        end if;
    end for;
    return rows;
end function;

verify_primes := [
    211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
    293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
    389,397,401,409,419,421,431,433,439,443
];

function ThreeBound(f)
    C := HyperellipticCurve(f); used := [];
    for p in verify_primes do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            nJ := Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            Append(~used,<p,nJ>);
            if nJ mod 3 ne 0 then return false,used; end if;
        catch e
            continue;
        end try;
    end for;
    return #used gt 0,used;
end function;

function SimpleWitness(f)
    C := HyperellipticCurve(f); R<T> := PolynomialRing(Q);
    for p in verify_primes do
      try
        fp := ChangeRing(f,GF(p));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
        lp := LPolynomial(HyperellipticCurve(fp));
        chi := T^4+(Q!Coefficient(lp,1))*T^3+(Q!Coefficient(lp,2))*T^2
                  +(Q!Coefficient(lp,3))*T+(Q!Coefficient(lp,4));
        if not IsIrreducible(chi) then continue; end if;
        K<pi> := NumberField(chi);
        if &and[Degree(MinimalPolynomial(pi^n)) eq 4:n in [2..12]] then
            return true,p,chi;
        end if;
      catch e
        continue;
      end try;
    end for;
    return false,0,R!0;
end function;

function ContactFromOrderThree(D,f)
    seq := Eltseq(D); q := seq[1]; vv := seq[2];
    if Degree(q) ne 2 then return false,q,vv,P!0; end if;
    A<m,n> := PolynomialRing(Q,2); PA<X> := PolynomialRing(A);
    phi := hom<P->PA|X>;
    qA := phi(q); vA := phi(vv); fA := phi(f);
    hA := vA+qA*(m*X+n);
    delta := hA^2-fA-m^2*qA^3;
    eqs := [Coefficient(delta,i):i in [0..Degree(delta)]];
    pts := Variety(ideal<A|eqs>);
    if #pts eq 0 then return false,q,vv,P!0; end if;
    mm := Q!pts[1][1]; nn := Q!pts[1][2];
    h := vv+q*(mm*x+nn);
    assert h^2-f eq mm^2*q^3;
    return true,q,vv,h;
end function;

rows := ReadRows(input_file); reduction_survivors := 0; exact := 0; hits := 0;
print "A2228_P1_PLUS3_VERIFY_START","file",input_file,"rows",#rows;
for rec in rows do
    family,n,d := Explode(rec); t := Q!n/d;
    vals := [];
    try vals := FamilyTuple(family,t); catch e continue; end try;
    if not Usable(vals) then continue; end if;
    f := CurvePolynomial(vals);
    survives,used := ThreeBound(f);
    if not survives then continue; end if;
    reduction_survivors +:= 1;
    print "EXTENDED_REDUCTION_SURVIVOR",rec,"used",used;
    if exact ge max_exact then continue; end if;
    exact +:= 1;
    C := HyperellipticCurve(f); J := Jacobian(C);
    G,mp := TorsionSubgroup(J); inv := Invariants(G);
    print "EXACT",rec,"t",t,"torsion",inv,"curve",f;
    if #G mod 3 ne 0 then continue; end if;
    hits +:= 1; D3 := J!0; found3 := false;
    for g in OrderedGenerators(G) do
        DD := mp(g);
        ord := Order(g);
        if ord mod 3 eq 0 then
            D3 := (ord div 3)*DD; found3 := true; break;
        end if;
    end for;
    contact := false; q := P!0; vv := P!0; h := P!0;
    if found3 then contact,q,vv,h := ContactFromOrderThree(D3,f); end if;
    simple,pcert,chi := SimpleWitness(f);
    print "TARGET_22224_HIT","record",rec,"t",t,"vals",vals,
          "torsion",inv,"D3",D3,"contact",contact,
          "q",q,"mumford_v",vv,"h",h,
          "simple",simple,"prime",pcert,"chi",chi;
end for;
print "A2228_P1_PLUS3_VERIFY_DONE","rows",#rows,
      "reduction_survivors",reduction_survivors,"exact",exact,"hits",hits;
UnsetLogFile();
quit;
