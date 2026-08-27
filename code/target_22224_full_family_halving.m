//////////////////////////////////////////////////////////////////////
// Full-family attack on [2,2,2,24].
//
// The marked direct A(2,2,2,12) construction has
//
//     D12 = D3 + D4,   ord(D3)=3,   ord(D4)=4.
//
// Since D3=2*(2D3), D12 is divisible by 2 if and only if D4 is.
// Thus the full (not point-on-the-curve) halving problem is the
// intersection A(2,2,2,8) + rational 3-torsion.  For the signed square
// branch model
//
//   y^2=x(x+a^2)(x+b^2)(x+c^2)(x+d^2),
//
// the marked A(2,2,2,8) cover is given by the four square conditions
//
//   abcd,
//   a(a+b)(a+c)(a+d),
//   b(b+a)(b+c)(b+d),
//   c(c+a)(c+b)(c+d).
//
// This is the full three-dimensional finite cover, not the
// two-dimensional M(2,2,2,8) K3 divisor s2^2=4s4 where a half happens
// to be represented by a point of the curve.
//
// Modes:
//   finite -- exhaustive finite-field target masks for this full cover;
//   bank   -- scan the entire supplied tor2228.txt integer bank for 3-primary
//             torsion, with exact verification of every survivor;
//   deep   -- exact audit of the two primitive bank points surviving the
//             contact-open mod-13 and forced mod-169 boundary conditions;
//   all    -- both.
//
// Typical runs from torsion_jac:
//
//   magma -b mode:=finite PrimeList:=11,13,17,19,23 \
//     code/target_22224_full_family_halving.m
//
//   magma -b mode:=bank \
//     code/target_22224_full_family_halving.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned mode then mode := "all"; end if;
if not assigned PrimeList then
    PrimeList := [7,11,13,17,19,23,29,31];
elif Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s) : s in Split(PrimeList,",")];
end if;
if not assigned bank_file then
    bank_file := "data/tor2228_bank.txt";
end if;
if not assigned max_rows then max_rows := 0;
elif Type(max_rows) eq MonStgElt then max_rows := StringToInteger(max_rows);
end if;
if not assigned log_file then
    log_file := "results/target_22224_full_family_halving.log";
end if;
if not assigned finite_output_prefix then
    finite_output_prefix := "results/target_22224_full_family_halving";
end if;

SetLogFile(log_file : Overwrite := true);

Q := Rationals();
Z := Integers();
PX<X> := PolynomialRing(Q);

function ReadTupleFile(filename)
    rows := Split(Read(filename),"\n");
    out := [];
    for raw in rows do
        row := raw;
        if #row gt 0 and row[#row] eq "\r" then row := row[1..#row-1]; end if;
        if #row lt 2 or row[1] ne "[" or row[#row] ne "]" then continue; end if;
        parts := Split(row[2..#row-1],",");
        if #parts eq 4 then Append(~out,[StringToInteger(s):s in parts]); end if;
    end for;
    return out;
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

function FullCoverQ(vals)
    return &and[r ne 0 and IsSquare(Q!r) : r in CoverRadicands(vals)];
end function;

function CurveFromTuple(vals)
    return X*&*[X+(Q!z)^2:z in vals];
end function;

// Coefficient keys of all smooth degree-two rational 3-torsion contact
// presentations over F_p.  This is the eliminated form of
// h^2-f=m^2 q^3, with V=v^2 and L=1/m.
function ContactCoefficientKeys(p)
    F := GF(p); PF<x> := PolynomialRing(F);
    keys := {}; witnesses := AssociativeArray(); raw := 0;
    for L in F do
      if L eq 0 then continue; end if;
      M := L^2;
      for U in F do for vv in F do
        // Do not discard vv=0, double-root q, or q meeting f.  These are
        // legitimate reductions of rational 3-torsion classes (the record
        // itself lies on the double-root q chart).
        q := x^2+U*x+vv^2;
        for e1 in F do
            // If P=8*L*Coefficient(h,1), coefficient comparison gives
            // P=4*M*e1+12*(U^2+vv^2)-(M+3*U)^2 and AA=P/2.
            // The original version accidentally doubled the final square.
            AA := (4*M*e1+12*(U^2+vv^2)-(M+3*U)^2)/2;
            e2 := ((M+3*U)*AA+8*vv^3-4*U^3-24*U*vv^2)/(4*M);
            e3 := (AA^2+16*(M+3*U)*vv^3
                    -48*(U^2*vv^2+vv^4))/(16*M);
            e4 := (AA*vv^3-6*U*vv^4)/(2*M);
            if e4 eq 0 then continue; end if;
            f := x^5+e1*x^4+e2*x^3+e3*x^2+e4*x;
            if Discriminant(f) eq 0 then continue; end if;
            key := <Z!e1,Z!e2,Z!e3,Z!e4>;
            Include(~keys,key); raw +:= 1;
            if not IsDefined(witnesses,key) then
                witnesses[key] := <Z!L,Z!U,Z!vv>;
            end if;
        end for;
      end for; end for;
    end for;
    return keys,witnesses,raw;
end function;

function IsProjectiveRepresentative(vals)
    first := 0;
    for i in [1..4] do if vals[i] ne 0 then first := i; break; end if; end for;
    return first ne 0 and vals[first] eq 1;
end function;

procedure FiniteMask(p)
    F := GF(p);
    keys,witnesses,raw := ContactCoefficientKeys(p);
    projective := 0; cover := 0; open_cover := 0; target := 0;
    curvekeys := {}; targetkeys := {}; samples := [];
    rank_hist := [0:i in [0..7]];
    // Jacobian of the full four-cover-plus-three-contact incidence.
    RR<aa0,bb0,cc0,dd0,r00,r10,r20,r30,LL0,UU0,vv0> :=
        PolynomialRing(F,11);
    base0 := [aa0,bb0,cc0,dd0];
    rad0 := CoverRadicands(base0);
    A0,B0,C0,D0 := Explode([z^2:z in base0]);
    ee1 := A0+B0+C0+D0;
    ee2 := A0*B0+A0*C0+A0*D0+B0*C0+B0*D0+C0*D0;
    ee3 := A0*B0*C0+A0*B0*D0+A0*C0*D0+B0*C0*D0;
    ee4 := A0*B0*C0*D0;
    MM0 := LL0^2;
    aux0 := (4*MM0*ee1+12*(UU0^2+vv0^2)
             -(MM0+3*UU0)^2)/2;
    incidence := [r00^2-rad0[1],r10^2-rad0[2],
                  r20^2-rad0[3],r30^2-rad0[4],
                  (MM0+3*UU0)*aux0+8*vv0^3-4*ee2*MM0
                    -4*UU0^3-24*UU0*vv0^2,
                  aux0^2+16*(MM0+3*UU0)*vv0^3-16*ee3*MM0
                    -48*(UU0^2*vv0^2+vv0^4),
                  aux0*vv0^3-2*ee4*MM0-6*UU0*vv0^4];
    Jac := Matrix(RR,7,11,[Derivative(e,j):e in incidence,j in [1..11]]);
    filename := Sprintf("%o_p%o.tsv",finite_output_prefix,p);
    out := Open(filename,"w");
    fprintf out,"a\tb\tc\td\te1\te2\te3\te4\tL\tU\tv\n";
    for a in F do for b in F do for c in F do for d in F do
        vals := [a,b,c,d];
        if not IsProjectiveRepresentative(vals) then continue; end if;
        projective +:= 1;
        rad := CoverRadicands(vals);
        if not &and[r ne 0 and IsSquare(r):r in rad] then continue; end if;
        cover +:= 1;
        sq := [z^2:z in vals];
        if &or[z eq 0:z in vals] or #Set(sq) ne 4 then continue; end if;
        open_cover +:= 1;
        A,B,C,D := Explode(sq);
        e1 := A+B+C+D;
        e2 := A*B+A*C+A*D+B*C+B*D+C*D;
        e3 := A*B*C+A*B*D+A*C*D+B*C*D;
        e4 := A*B*C*D;
        key := <Z!e1,Z!e2,Z!e3,Z!e4>;
        Include(~curvekeys,key);
        if key notin keys then continue; end if;
        target +:= 1; Include(~targetkeys,key);
        W := witnesses[key];
        roots := [SquareRoot(z):z in rad];
        coords := vals cat roots cat [F!W[1],F!W[2],F!W[3]];
        jr := Rank(Matrix(F,7,11,[Evaluate(Jac[i,j],coords)
                                  :i in [1..7],j in [1..11]]));
        rank_hist[jr+1] +:= 1;
        fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            Z!a,Z!b,Z!c,Z!d,Z!e1,Z!e2,Z!e3,Z!e4,W[1],W[2],W[3];
        if #samples lt 8 then Append(~samples,<vals,key,W>); end if;
    end for; end for; end for; end for;
    delete out;
    print "FULL_FAMILY_FINITE","p",p,"contact_keys",#keys,
          "raw_contact_parameters",raw,"projective_bases",projective,
          "cover_bases",cover,"open_cover_bases",open_cover,
          "open_curvekeys",#curvekeys,"target_bases",target,
          "target_curvekeys",#targetkeys,"incidence_rank_hist_0_to_7",
          rank_hist,"mask_file",filename;
    print "FULL_FAMILY_FINITE_SAMPLES",p,samples;
end procedure;

bank_primes := [
    5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,
    149,151,157,163,167,173,179,181,191,193,197,199
];

function ThreePrimaryBound(f)
    C := HyperellipticCurve(f); used := []; g := 0;
    for p in bank_primes do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            n := Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g := (g eq 0) select n else GCD(g,n);
            Append(~used,<p,n,g>);
            if g mod 3 ne 0 then return false,g,used; end if;
        catch e
            continue;
        end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

function RootPowerWitness(f)
    C := HyperellipticCurve(f);
    for p in bank_primes do
      try
        fp := ChangeRing(f,GF(p));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
        lp := LPolynomial(HyperellipticCurve(fp));
        R<T> := PolynomialRing(Q);
        chi := T^4+(Q!Coefficient(lp,1))*T^3
                    +(Q!Coefficient(lp,2))*T^2
                    +(Q!Coefficient(lp,3))*T+(Q!Coefficient(lp,4));
        if not IsIrreducible(chi) then continue; end if;
        K<pi> := NumberField(chi);
        degs := [Degree(MinimalPolynomial(pi^n)):n in [2..12]];
        if &and[z eq 4:z in degs] then return true,p,chi,degs; end if;
      catch e
        continue;
      end try;
    end for;
    R<T> := PolynomialRing(Q);
    return false,0,R!0,[];
end function;

procedure BankScan()
    rows := ReadTupleFile(bank_file);
    stop := (max_rows eq 0) select #rows else Min(#rows,max_rows);
    eliminated := AssociativeArray();
    for p in bank_primes do eliminated[p] := 0; end for;
    cover_ok := 0; survivors := 0; exact := 0; hits := 0;
    print "FULL_FAMILY_BANK_START","file",bank_file,"rows",#rows,"stop",stop;
    for i in [1..stop] do
        vals := rows[i];
        if not FullCoverQ(vals) then
            print "FULL_FAMILY_BANK_COVER_FAILURE",i,vals,CoverRadicands(vals);
            continue;
        end if;
        cover_ok +:= 1;
        f := CurveFromTuple(vals);
        survives,g,used := ThreePrimaryBound(f);
        if not survives then
            eliminated[used[#used][1]] +:= 1;
            continue;
        end if;
        survivors +:= 1;
        print "FULL_FAMILY_BANK_SURVIVOR",i,vals,"gcd_bound",g,"used",used;
        exact +:= 1;
        G,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
        inv := Invariants(G);
        print "FULL_FAMILY_BANK_EXACT",i,vals,"invariants",inv,"order",#G;
        if #G mod 3 ne 0 then continue; end if;
        hits +:= 1;
        simple,pcert,chi,degs := RootPowerWitness(f);
        print "TARGET_22224_HIT","row",i,"tuple",vals,"curve",f,
              "torsion",inv,"simple",simple,"prime",pcert,
              "chi",chi,"degrees",degs;
    end for;
    print "FULL_FAMILY_BANK_ELIMINATIONS",
          [<p,eliminated[p]>:p in bank_primes|eliminated[p] gt 0];
    print "FULL_FAMILY_BANK_DONE","rows",stop,"cover_verified",cover_ok,
          "reduction_survivors",survivors,"exact_tests",exact,"hits",hits;
end procedure;

procedure DeepSeedAudit()
    seeds := [
        [121,1919,3211,4949],
        [1369,1711,2349,9971]
    ];
    print "FULL_FAMILY_DEEP13_SEED_AUDIT_START","count",#seeds;
    for vals in seeds do
        assert FullCoverQ(vals);
        f := CurveFromTuple(vals);
        survives,g,used := ThreePrimaryBound(f);
        G,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
        simple,pcert,chi,degs := RootPowerWitness(f);
        print "FULL_FAMILY_DEEP13_SEED","tuple",vals,
              "three_bound_survives",survives,"gcd_bound",g,
              "first_killing_prime",used[#used][1],
              "torsion",Invariants(G),"order",#G,
              "strict_simple",simple,"certificate_prime",pcert,
              "chi",chi,"degrees",degs;
    end for;
    print "FULL_FAMILY_DEEP13_SEED_AUDIT_DONE";
end procedure;

print "TARGET_22224_FULL_FAMILY_HALVING_START","mode",mode,
      "PrimeList",PrimeList;
if mode in {"finite","all"} then for p in PrimeList do FiniteMask(p); end for; end if;
if mode in {"bank","all"} then BankScan(); end if;
if mode in {"deep","all"} then DeepSeedAudit(); end if;
print "TARGET_22224_FULL_FAMILY_HALVING_DONE";
UnsetLogFile();
quit;
