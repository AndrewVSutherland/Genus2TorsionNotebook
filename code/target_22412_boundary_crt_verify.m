//////////////////////////////////////////////////////////////////////
// Exact verifier and finite positive controls for the boundary/CRT
// attack on [2,2,4,12].
//
// The Python driver writes any smooth q-square-slice survivor as twelve
// rationals followed by its second-half partition labels.  This file
// reconstructs the direct A(2,2,2,12) contact identity, computes the exact
// torsion subgroup, and requires a root-power simplicity certificate for
// every target hit.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned candidate_file then
    candidate_file :=
      "results/target_22412_boundary_crt_candidates.txt";
end if;
if not assigned log_file then
    log_file :=
      "results/target_22412_boundary_crt_verify.log";
end if;

SetLogFile(log_file : Overwrite := true);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function RationalFromString(s)
    parts := Split(s,"/");
    assert #parts eq 2;
    return Q!StringToInteger(parts[1])/StringToInteger(parts[2]);
end function;

function DirectData(K,aa,bb,cc,dd,uu,tt,vv)
    PK<X> := PolynomialRing(K);
    f := X*(X+aa^2)*(X+bb^2)*(X+cc^2)*(X+dd^2);
    q := X^2+uu*X+tt^2;
    h := X^3+(1+3*uu)/2*X^2+vv*X+tt^3;
    return PK,f,q,h;
end function;

function SecondHalfRadicands(squares,ord)
    AA:=squares[ord[1]]; BB:=squares[ord[2]];
    CC:=squares[ord[3]]; DD:=squares[ord[4]];
    return [
      (CC-AA)*(CC-BB),
      (CC-AA)*(DD-AA),
      (CC-BB)*(DD-BB),
      (DD-AA)*(DD-BB)
    ];
end function;

function IsDivisibleBy2Finite(J,D)
    GG,phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(GG);
    for i in [1..#coords] do
        if GCD(2,invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function IntegralSquareModel(f)
    den := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return Parent(f)!(den^2*f),den;
end function;

function FrobeniusPolynomial(C,p)
    Lp := LPolynomial(ChangeRing(C,GF(p)));
    R<T> := PolynomialRing(Q);
    return &+[Q!Coefficient(Lp,i)*T^(4-i) : i in [0..4]];
end function;

function RootPowerCertificate(f)
    C := HyperellipticCurve(f);
    for pp in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,
               59,61,67,71,73,79,83,89,97,101,103,107,109] do
        try
            fp := ChangeRing(f,GF(pp));
            if Degree(fp) notin {5,6} or Discriminant(fp) eq 0 then
                continue;
            end if;
            chi := FrobeniusPolynomial(C,pp);
            if not IsIrreducible(chi) then continue; end if;
            K<pi> := NumberField(chi);
            degrees := [Degree(MinimalPolynomial(pi^n)) : n in [2..12]];
            if &and[d eq 4 : d in degrees] then
                return true,pp,chi,degrees;
            end if;
        catch e
            continue;
        end try;
    end for;
    R<T> := PolynomialRing(Q);
    return false,0,R!0,[];
end function;

procedure GeneralFiniteControl(pp,vals,ord)
    K := GF(pp);
    PK<X> := PolynomialRing(K);
    roots := [K!z : z in vals];
    squares := [z^2 : z in roots];
    f := X*&*[X+s : s in squares];
    assert Discriminant(f) ne 0;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    rad := SecondHalfRadicands(squares,ord);
    assert &and[IsSquare(r) : r in rad];
    T2 := J![(X+squares[ord[1]])*(X+squares[ord[2]]),K!0];
    divisible := IsDivisibleBy2Finite(J,T2);
    has3 := #J mod 3 eq 0;
    print "GENERAL_POSITIVE_CONTROL",pp,"roots",roots,"partition",ord,
          "radicands",rad,"divisible",divisible,"J_order",#J,
          "three_primary",has3;
    assert divisible and has3;
end procedure;

procedure QSquareFiniteControl(pp,A0,B0,expected_d2,ord)
    K := GF(pp);
    AA := K!A0; BB := K!B0; CC := 1/(AA*BB);
    RR := AA^2+BB^2+CC^2-3;
    SS := 1/AA^2+1/BB^2+1/CC^2-3;
    assert RR ne 0 and SS ne 0 and IsSquare(RR) and IsSquare(SS);
    rho := SquareRoot(RR); sigma := SquareRoot(SS);
    DD2 := SS/RR;
    assert DD2 eq K!expected_d2 and IsSquare(DD2);
    DD := SquareRoot(DD2);

    scale := 1/(2*rho); tt := scale^2; uu := 2*tt;
    aa := scale*AA; bb := scale*BB; cc := scale*CC;
    dd := 2*scale^2*sigma;
    vv := 3*tt^2+dd^2/2;
    PK,f,q,h := DirectData(K,aa,bb,cc,dd,uu,tt,vv);
    assert h^2-f eq q^3 and Discriminant(f) ne 0;
    J := Jacobian(HyperellipticCurve(f));
    D3 := J![q,h mod q];
    squares := [AA^2,BB^2,CC^2,DD2];
    rad := SecondHalfRadicands(squares,ord);
    assert &and[IsSquare(r) : r in rad];
    actual_squares := [aa^2,bb^2,cc^2,dd^2];
    T2 := J![(PK.1+actual_squares[ord[1]])*
             (PK.1+actual_squares[ord[2]]),K!0];
    divisible := IsDivisibleBy2Finite(J,T2);
    print "QSQUARE_POSITIVE_CONTROL",pp,"A",AA,"B",BB,"C",CC,
          "D2",DD2,"partition",ord,"D3_order",Order(D3),
          "second_half",divisible,"J_order",#J;
    assert Order(D3) eq 3 and divisible and #J mod 3 eq 0;
end procedure;

print "TARGET_22412_BOUNDARY_CRT_VERIFY_START";

// Positive controls from the global standard cover at p=17,19.
GeneralFiniteControl(17,[1,4,2,8],[1,2,3,4]);
GeneralFiniteControl(19,[1,4,3,7],[1,2,3,4]);

// The q-square slice itself first has smooth second-half controls at p=43.
QSquareFiniteControl(43,2,3,15,[1,4,2,3]);
QSquareFiniteControl(47,2,17,37,[1,2,3,4]);

raw := Read(candidate_file);
rows := Split(raw,"\n");
checked := 0;
target_hits := 0;
for row in rows do
    if #row eq 0 or row[1] eq "#" then continue; end if;
    fields := Split(row," ");
    assert #fields eq 13;
    z := [RationalFromString(fields[i]) : i in [1..12]];
    A0:=z[1]; B0:=z[2]; C0:=z[3]; rho:=z[4]; sigma:=z[5];
    aa:=z[6]; bb:=z[7]; cc:=z[8]; dd:=z[9];
    uu:=z[10]; tt:=z[11]; vv:=z[12];
    assert A0*B0*C0 eq 1;
    assert A0^2+B0^2+C0^2-3 eq rho^2;
    assert 1/A0^2+1/B0^2+1/C0^2-3 eq sigma^2;
    PK,f,q,h := DirectData(Q,aa,bb,cc,dd,uu,tt,vv);
    assert h^2-f eq q^3 and Discriminant(f) ne 0;

    fI,scaleI := IntegralSquareModel(P!f);
    TG,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
    invs := Invariants(TG);
    checked +:= 1;
    print "EXACT_CANDIDATE",checked,"A",A0,"B",B0,"C",C0,
          "partitions",fields[13],"torsion",invs,"order",#TG,
          "curve",fI;
    if invs eq [2,2,4,12] then
        target_hits +:= 1;
        cert,pcert,chi,degrees := RootPowerCertificate(fI);
        print "TARGET_22412_EXACT_HIT","simple",cert,"prime",pcert,
              "chi",chi,"degrees",degrees;
        assert cert;
    end if;
end for;

print "EXACT_CANDIDATE_SUMMARY","checked",checked,"target_hits",target_hits;
print "TARGET_22412_BOUNDARY_CRT_VERIFY_DONE";
UnsetLogFile();
quit;
