//////////////////////////////////////////////////////////////////////
//  Target [2,2,4,12] from A(2,2,2,12): a second independent half
//  of a rational J[2] class.
//
//  This driver has three parts.
//
//  (1) record:
//      Enumerate all 15 nonzero J[2](Q) classes of the order-96
//      [2,2,2,12] record.  Compare the exact Stoll--Zarhin
//      squareclass criterion with Magma's IsDivisibleBy.
//
//  (2) finite:
//      On the standard labelled A(2,2,4,4) cover
//
//        y^2=x(x+a^2)(x+b^2)(x+c^2)(x+d^2),
//
//      count projective residue classes, open second-half points, and
//      open points for which 3 divides #J(F_p).
//
//  (3) bank:
//      Scan the existing 26,653-row A(2,2,4,4) tuple bank for the
//      necessary rational-3 reduction condition.  Exact-test every
//      survivor.  Any target hit gets a full root-power simplicity
//      check, not merely an irreducible Frobenius quartic.
//
//  Typical complete run from the repository root:
//
//    magma -b mode:=all \
//      code/target_22412_second_j2_halving.m
//
//  Individual modes use mode:=record, mode:=finite, or mode:=bank.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned mode then
    mode := "all";
end if;
if not assigned input_file then
    input_file := "data/tor2244_bank.txt";
end if;
if not assigned max_rows then
    max_rows := 0;
elif Type(max_rows) eq MonStgElt then
    max_rows := StringToInteger(max_rows);
end if;
if not assigned max_exact then
    max_exact := 50;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 5000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned finite_prime_max then
    finite_prime_max := 19;
elif Type(finite_prime_max) eq MonStgElt then
    finite_prime_max := StringToInteger(finite_prime_max);
end if;
if not assigned validate_limit then
    validate_limit := 4;
elif Type(validate_limit) eq MonStgElt then
    validate_limit := StringToInteger(validate_limit);
end if;
if not assigned write_log then
    write_log := true;
elif Type(write_log) eq MonStgElt then
    write_log := write_log in {"true","True","1","yes"};
end if;
if not assigned log_file then
    log_file := "results/target_22412_second_j2_halving.log";
end if;

if write_log then
    SetLogFile(log_file : Overwrite := true);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

//////////////////////////////////////////////////////////////////////
// General exact squareclass criterion.
//////////////////////////////////////////////////////////////////////

function AllPairs6()
    return &cat[[<i,j> : j in [i+1..6]] : i in [1..5]];
end function;

function PairPolynomial(roots,ij)
    return (x-Q!roots[ij[1]])*(x-Q!roots[ij[2]]);
end function;

// For the class represented by the branch pair {A,B}, send A to
// infinity and B to zero.  If c runs over the other four roots, put
//
//   D = product_c(A-c),
//   q_c = -D*(c-B)/(c-A).
//
// The transformed odd model is Y^2=X*product_c(X+q_c), so the class
// is twice a rational class iff all four q_c are rational squares.
function ExactRadicands(roots,ij)
    A := Q!roots[ij[1]];
    B := Q!roots[ij[2]];
    rem := [k : k in [1..6] | k ne ij[1] and k ne ij[2]];
    cs := [Q!roots[k] : k in rem];
    DD := &*[A-c : c in cs];
    vals := [-DD*(c-B)/(c-A) : c in cs];
    return vals,rem,DD;
end function;

function AllSquaresQ(vals)
    return &and[v ne 0 and IsSquare(Q!v) : v in vals];
end function;

function RationalSquareRoots(vals)
    out := [];
    for v in vals do
        ok,s := IsSquare(Q!v);
        assert ok;
        Append(~out,s);
    end for;
    return out;
end function;

function PairIntersectionParity(ij,kl)
    return (#[a : a in [ij[1],ij[2]] | a in {kl[1],kl[2]}]) mod 2;
end function;

// Standard A(2,2,4,4) cover after the already-divisible pair has
// been moved to {0,infinity}.  The entries are A=a^2, ..., D=d^2.
// The ordered partition is [A,B | C,D].
function SecondHalfRadicands(squares,ord)
    AA := squares[ord[1]];
    BB := squares[ord[2]];
    CC := squares[ord[3]];
    DD := squares[ord[4]];
    return [
        (CC-AA)*(CC-BB),
        (CC-AA)*(DD-AA),
        (CC-BB)*(DD-BB),
        (DD-AA)*(DD-BB)
    ];
end function;

partition_orders := [
    [1,2,3,4],
    [1,3,2,4],
    [1,4,2,3]
];

function SquareclassRankQ(vals)
    primes := [];
    for vv in vals do
        v := Q!vv;
        assert v ne 0;
        for pe in Factorization(Abs(Numerator(v))) cat
                  Factorization(Denominator(v)) do
            if pe[1] notin primes then Append(~primes,pe[1]); end if;
        end for;
    end for;
    Sort(~primes);
    F2 := GF(2);
    M := ZeroMatrix(F2,#vals,#primes+1);
    for i in [1..#vals] do
        v := Q!vals[i];
        if v lt 0 then M[i,1] := 1; end if;
        for j in [1..#primes] do
            pp := primes[j];
            M[i,j+1] := F2!(Valuation(Numerator(v),pp)
                - Valuation(Denominator(v),pp));
        end for;
    end for;
    return Rank(M),M,primes;
end function;

//////////////////////////////////////////////////////////////////////
// Exact torsion and geometric-simplicity helpers.
//////////////////////////////////////////////////////////////////////

function ContainsTarget22412(invs)
    if #invs gt 4 then return false; end if;
    pad := [1 : i in [1..4-#invs]] cat [Z!n : n in invs];
    target := [2,2,4,12];
    return &and[pad[i] mod target[i] eq 0 : i in [1..4]];
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
            if Degree(fp) ne Degree(f) or Discriminant(fp) eq 0 then
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

//////////////////////////////////////////////////////////////////////
// Part 1: source record audit and symbolic global cover.
//////////////////////////////////////////////////////////////////////

procedure RunRecordAudit()
    print "TARGET_22412_RECORD_AUDIT_START";

    roots := [
        Q!(-75593)/145,
        Q!(-519),
        Q!(-39),
        Q!(-299)/12,
        Q!59,
        Q!123109/1740
    ];
    lc := Q!3027600;
    assert IsSquare(lc);
    f := lc*&*[x-r : r in roots];
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    G,mp := TorsionSubgroup(J);
    print "RECORD_ROOTS",roots;
    print "RECORD_TORSION",Invariants(G),"order",#G;
    assert Invariants(G) eq [2,2,2,12];

    pairs := AllPairs6();
    divisible_pairs := [];
    for k in [1..#pairs] do
        ij := pairs[k];
        vals,rem,DD := ExactRadicands(roots,ij);
        criterion := AllSquaresQ(vals);
        T2 := J![PairPolynomial(roots,ij),Q!0];
        divisible,half := IsDivisibleBy(T2,2);
        assert criterion eq divisible;
        print "RECORD_J2_CLASS",k,"pair",ij,
              "roots",[roots[ij[1]],roots[ij[2]]],
              "u",PairPolynomial(roots,ij),
              "radicands",vals,
              "criterion",criterion,"magma",divisible;
        if divisible then
            Append(~divisible_pairs,ij);
            print "RECORD_J2_HALF","pair",ij,"half_order",Order(half),
                  "check",2*half eq T2,
                  "square_roots",RationalSquareRoots(vals);
        end if;
    end for;
    print "RECORD_DIVISIBLE_PAIRS",divisible_pairs;
    assert #divisible_pairs eq 1;
    known := divisible_pairs[1];
    assert PairPolynomial(roots,known) eq (x+519)*(x+39);

    cross := [];
    disjoint := [];
    for ij in pairs do
        if ij eq known then continue; end if;
        if PairIntersectionParity(known,ij) eq 1 then
            Append(~cross,ij);
        else
            Append(~disjoint,ij);
        end if;
    end for;
    print "LEVEL2_STABILIZER_ORBITS",
          "known",known,"cross_nonisotropic_count",#cross,
          "cross",cross,"disjoint_isotropic_count",#disjoint,
          "disjoint",disjoint;
    assert #cross eq 8 and #disjoint eq 6;

    // Move the known pair A,B to infinity,zero.  Its four exact
    // radicands are squares and become the A,B,C,D in the standard
    // odd model x*product(x+A_i).
    std_squares,rem,DD := ExactRadicands(roots,known);
    std_roots := RationalSquareRoots(std_squares);
    std_normalized := [u/std_roots[1] : u in std_roots];
    print "RECORD_STANDARD_ODD_SQUARES",std_squares;
    print "RECORD_STANDARD_ODD_SQUARE_ROOTS",std_roots;
    print "RECORD_STANDARD_ODD_NORMALIZED_ROOTS",std_normalized;
    print "RECORD_STANDARD_ODD_POLYNOMIAL",x*&*[x+s : s in std_squares];
    assert AllSquaresQ(std_squares);
    assert std_normalized eq [
        Q!1,Q!169/14399,Q!91/3179,Q!5681/189431
    ];

    second_partitions := [];
    for ord in partition_orders do
        vals := SecondHalfRadicands(std_squares,ord);
        ok := AllSquaresQ(vals);
        print "RECORD_SECOND_HALF_PARTITION",ord,
              "radicands",vals,"all_squares",ok;
        if ok then Append(~second_partitions,ord); end if;
    end for;
    print "RECORD_SECOND_HALF_PARTITIONS",second_partitions;
    assert #second_partitions eq 0;

    // Symbolic form of the single level-2 orbit of possible second
    // independent halves.  The product relation leaves three generic
    // squareclasses, hence a degree-8 labelled cover.
    SR<AA,BB,CC,DDs> := PolynomialRing(Q,4);
    symbolic := [
        (CC-AA)*(CC-BB),
        (CC-AA)*(DDs-AA),
        (CC-BB)*(DDs-BB),
        (DDs-AA)*(DDs-BB)
    ];
    assert &*symbolic eq
        ((CC-AA)*(CC-BB)*(DDs-AA)*(DDs-BB))^2;
    sample_vals := [Q!Evaluate(g,[1,4,9,16]) : g in symbolic];
    scrank,M,primes := SquareclassRankQ(sample_vals);
    assert scrank eq 3;
    print "GLOBAL_SECOND_HALF_COVER_RADICANDS",symbolic;
    print "GLOBAL_SECOND_HALF_PRODUCT_RELATION",
          "product = ((C-A)(C-B)(D-A)(D-B))^2";
    print "GLOBAL_SECOND_HALF_SAMPLE",[1,4,9,16],sample_vals,
          "squareclass_rank",scrank,"generic_degree",2^scrank;
    print "GLOBAL_A22212_COUPLING",
          "h^2-q^3=x*(x+A)*(x+B)*(x+C)*(x+D),",
          "A=a^2,B=b^2,C=c^2,D=d^2, plus the four radicand squares";
    print "TARGET_22412_RECORD_AUDIT_DONE";
end procedure;

//////////////////////////////////////////////////////////////////////
// Part 2: finite-field residue/component analysis.
//////////////////////////////////////////////////////////////////////

function IsProjectiveRepresentative(vals)
    first := 0;
    for i in [1..#vals] do
        if vals[i] ne 0 then first := i; break; end if;
    end for;
    if first eq 0 then return false; end if;
    return vals[first] eq 1;
end function;

function OpenSquareTuple(vals)
    if &or[a eq 0 : a in vals] then return false; end if;
    sq := [a^2 : a in vals];
    return #Setseq(Seqset(sq)) eq 4;
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

procedure FinitePrimeAnalysis(pp)
    F := GF(pp);
    PF<X> := PolynomialRing(F);
    projective := 0;
    cover_total := 0;
    cover_boundary := 0;
    open_total := 0;
    open_cover := 0;
    open_cover_3 := 0;
    validated := 0;
    validated_positive := 0;
    validated_negative := 0;
    boundary_signatures := AssociativeArray();
    first_open_cover := [];
    first_open_cover_3 := [];

    for a in F do
      for b in F do
       for c in F do
        for d in F do
            vals := [a,b,c,d];
            if not IsProjectiveRepresentative(vals) then continue; end if;
            projective +:= 1;
            sq := [z^2 : z in vals];
            rad := SecondHalfRadicands(sq,[1,2,3,4]);
            cover := &and[IsSquare(r) : r in rad];
            open := OpenSquareTuple(vals);
            if open then open_total +:= 1; end if;
            if cover then
                cover_total +:= 1;
                if open then
                    open_cover +:= 1;
                    if #first_open_cover lt 3 then Append(~first_open_cover,vals); end if;
                else
                    cover_boundary +:= 1;
                    zeros := #[z : z in vals | z eq 0];
                    collisions := 0;
                    for i in [1..3] do
                        for j in [i+1..4] do
                            if sq[i] eq sq[j] then collisions +:= 1; end if;
                        end for;
                    end for;
                    key := Sprintf("zeros_%o_collisions_%o",zeros,collisions);
                    if IsDefined(boundary_signatures,key) then
                        boundary_signatures[key] +:= 1;
                    else
                        boundary_signatures[key] := 1;
                    end if;
                end if;
            end if;

            if not open then continue; end if;
            f := X*&*[X+s : s in sq];
            assert Degree(f) eq 5 and Discriminant(f) ne 0;
            C := HyperellipticCurve(f);
            J := Jacobian(C);

            // Validate the first few open tuples, but also force a
            // positive and a negative test whenever both kinds occur.
            // This prevents a small validate_limit from checking only
            // the lexicographically earliest (usually negative) cases.
            need_validation := validated lt validate_limit or
                (cover and validated_positive eq 0) or
                ((not cover) and validated_negative eq 0);
            if need_validation then
                T2 := J![(X+sq[1])*(X+sq[2]),F!0];
                assert IsDivisibleBy2Finite(J,T2) eq cover;
                validated +:= 1;
                if cover then
                    validated_positive +:= 1;
                else
                    validated_negative +:= 1;
                end if;
            end if;

            if cover and (#J mod 3 eq 0) then
                open_cover_3 +:= 1;
                if #first_open_cover_3 lt 3 then
                    Append(~first_open_cover_3,vals);
                end if;
            end if;
        end for;
       end for;
      end for;
    end for;

    print "FINITE_PRIME",pp,
          "projective",projective,
          "open_total",open_total,
          "cover_total",cover_total,
          "cover_boundary",cover_boundary,
          "open_cover",open_cover,
          "open_cover_with_3",open_cover_3,
          "validated",validated,
          "validated_positive",validated_positive,
          "validated_negative",validated_negative;
    print "FINITE_OPEN_COVER_SAMPLES",pp,first_open_cover;
    print "FINITE_OPEN_COVER_3_SAMPLES",pp,first_open_cover_3;
    print "FINITE_BOUNDARY_SIGNATURES",pp,
          Sort([<k,boundary_signatures[k]> : k in Keys(boundary_signatures)]);
end procedure;

procedure RunFiniteAnalysis()
    print "TARGET_22412_FINITE_ANALYSIS_START";
    for pp in [p : p in [5..finite_prime_max] | IsPrime(p) and p ne 2 and p ne 3] do
        FinitePrimeAnalysis(pp);
    end for;
    print "TARGET_22412_FINITE_ANALYSIS_DONE";
end procedure;

//////////////////////////////////////////////////////////////////////
// Part 3: existing A(2,2,4,4) rational tuple bank plus 3-torsion.
//////////////////////////////////////////////////////////////////////

function ReadTupleFile(filename)
    S := Read(filename);
    rows := Split(S,"\n");
    tuples := [];
    for raw in rows do
        row := raw;
        if #row gt 0 and row[#row] eq "\r" then
            row := row[1..#row-1];
        end if;
        if #row ge 2 and row[1] eq "[" and row[#row] eq "]" then
            parts := Split(row[2..#row-1],",");
            tup := [StringToInteger(z) : z in parts];
            if #tup eq 4 then Append(~tuples,tup); end if;
        end if;
    end for;
    return tuples;
end function;

function TuplePolynomial(tup)
    return x*&*[x+(Q!a)^2 : a in tup];
end function;

function TupleCoverPartitions(tup)
    sq := [(Q!a)^2 : a in tup];
    out := [];
    for ord in partition_orders do
        vals := SecondHalfRadicands(sq,ord);
        if AllSquaresQ(vals) then Append(~out,ord); end if;
    end for;
    return out;
end function;

function PassesRationalThreeReductions(f,primes)
    used := [];
    for pp in primes do
        fp := ChangeRing(f,GF(pp));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
        n := Z!#Jacobian(HyperellipticCurve(fp));
        Append(~used,<pp,n>);
        if n mod 3 ne 0 then return false,pp,used; end if;
    end for;
    return true,0,used;
end function;

// Boundary components in the labelled square-root coordinates.  Zi is
// a_i=0; Eij+ and Eij- are a_i=a_j and a_i=-a_j respectively.
function TupleBoundaryLabelsModp(tup,pp)
    F := GF(pp);
    vals := [F!a : a in tup];
    labels := [];
    for i in [1..4] do
        if vals[i] eq 0 then Append(~labels,Sprintf("Z%o",i)); end if;
    end for;
    for i in [1..3] do
        for j in [i+1..4] do
            if vals[i] eq vals[j] then
                Append(~labels,Sprintf("E%o%o+",i,j));
            end if;
            if vals[i] eq -vals[j] then
                Append(~labels,Sprintf("E%o%o-",i,j));
            end if;
        end for;
    end for;
    return labels;
end function;

procedure IncrementLabelCounts(~counts,labels)
    for label in labels do
        if IsDefined(counts,label) then
            counts[label] +:= 1;
        else
            counts[label] := 1;
        end if;
    end for;
end procedure;

procedure RunBankSearch()
    print "TARGET_22412_BANK_SEARCH_START";
    tuples := ReadTupleFile(input_file);
    if max_rows gt 0 and max_rows lt #tuples then
        tuples := tuples[1..max_rows];
    end if;
    primes := [5,7,11,13,17,19,23,29,31,37,41,43,47,53,
               59,61,67,71,73,79,83,89,97];
    print "BANK_CONFIG","input",input_file,"rows",#tuples,
          "primes",primes,"max_exact",max_exact;

    checked := 0;
    cover_verified := 0;
    reduction_survivors := 0;
    exact := 0;
    target_hits := 0;
    no_cover := 0;
    kill_counts := AssociativeArray(Z);
    boundary11 := AssociativeArray();
    boundary13 := AssociativeArray();
    good11 := 0;
    good13 := 0;
    for pp in primes do kill_counts[pp] := 0; end for;

    for tup in tuples do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "BANK_PROGRESS",checked,"cover_verified",cover_verified,
                  "reduction_survivors",reduction_survivors,
                  "exact",exact,"target_hits",target_hits;
        end if;
        parts := TupleCoverPartitions(tup);
        if #parts eq 0 then
            no_cover +:= 1;
            continue;
        end if;
        cover_verified +:= 1;
        labels11 := TupleBoundaryLabelsModp(tup,11);
        labels13 := TupleBoundaryLabelsModp(tup,13);
        if #labels11 eq 0 then
            good11 +:= 1;
        else
            IncrementLabelCounts(~boundary11,labels11);
        end if;
        if #labels13 eq 0 then
            good13 +:= 1;
        else
            IncrementLabelCounts(~boundary13,labels13);
        end if;
        f := TuplePolynomial(tup);
        assert Discriminant(f) ne 0;
        ok,badp,used := PassesRationalThreeReductions(f,primes);
        if not ok then
            kill_counts[badp] +:= 1;
            continue;
        end if;
        reduction_survivors +:= 1;
        print "BANK_REDUCTION_SURVIVOR",tup,
              "cover_partitions",parts,"used",used;
        if exact ge max_exact then continue; end if;

        C := HyperellipticCurve(f);
        TG,mp := TorsionSubgroup(Jacobian(C));
        invs := Invariants(TG);
        exact +:= 1;
        print "BANK_EXACT_TORSION",tup,invs,"order",#TG;
        if ContainsTarget22412(invs) then
            target_hits +:= 1;
            cert,pp,chi,degrees := RootPowerCertificate(f);
            print "TARGET_22412_HIT",tup,"torsion",invs,
                  "curve",f,"root_power_certificate",cert,
                  "prime",pp,"chi",chi,"degrees",degrees;
            assert invs ne [2,2,4,12] or cert;
        end if;
    end for;

    print "TARGET_22412_BANK_SEARCH_DONE",
          "checked",checked,"cover_verified",cover_verified,
          "no_cover",no_cover,
          "reduction_survivors",reduction_survivors,
          "exact",exact,"target_hits",target_hits;
    print "BANK_FIRST_KILL_COUNTS",
          [<pp,kill_counts[pp]> : pp in primes];
    print "BANK_BOUNDARY_COMPONENT_MARGINALS_P11","good",good11,
          Sort([<k,boundary11[k]> : k in Keys(boundary11)]);
    print "BANK_BOUNDARY_COMPONENT_MARGINALS_P13","good",good13,
          Sort([<k,boundary13[k]> : k in Keys(boundary13)]);
end procedure;

//////////////////////////////////////////////////////////////////////
// Dispatch.
//////////////////////////////////////////////////////////////////////

print "TARGET_22412_SECOND_J2_HALVING";
print "CONFIG","mode",mode,"finite_prime_max",finite_prime_max,
      "max_rows",max_rows,"max_exact",max_exact,
      "write_log",write_log,"log_file",log_file;

if mode eq "record" or mode eq "all" then RunRecordAudit(); end if;
if mode eq "finite" or mode eq "all" then RunFiniteAnalysis(); end if;
if mode eq "bank" or mode eq "all" then RunBankSearch(); end if;

print "TARGET_22412_SECOND_J2_HALVING_DONE";
quit;
