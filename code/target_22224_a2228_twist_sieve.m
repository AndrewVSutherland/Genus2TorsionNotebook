//////////////////////////////////////////////////////////////////////
// Quadratic-twist sieve for every non-bielliptic A(2,2,2,8) bank point.
//
// For z=(a,b,c,d), the curve is y^2=x*Prod(x+z_i^2).  A quadratic
// twist retains full rational 2-torsion and geometric simplicity, but may
// redistribute higher 2-power and odd torsion.  At each good prime there
// are only two local twist classes.  If neither Jacobian order is divisible
// by 192, no quadratic twist of that curve can realize [2,2,2,24].
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetLogFile("results/target_22224_a2228_twist_sieve.log" : Overwrite := true);
Z := Integers(); Q := Rationals(); P<x> := PolynomialRing(Q);

function Primitive(v)
    g := GCD([Z!Abs(z) : z in v]);
    if g gt 1 then v := [z div g : z in v]; end if;
    for z in v do
        if z ne 0 then
            if z lt 0 then v := [-w : w in v]; end if;
            break;
        end if;
    end for;
    return v;
end function;

function Rectangle(v)
    a,b,c,d := Explode([Abs(z) : z in v]);
    return a*b eq c*d or a*c eq b*d or a*d eq b*c;
end function;

function ReadBank(path)
    ans := {}; 
    for line in Split(Read(path),"\n") do
        if #line lt 3 or line[1] ne "[" or line[#line] ne "]" then continue; end if;
        ss := Split(line[2..#line-1],",");
        if #ss ne 4 then continue; end if;
        v := Primitive([StringToInteger(s) : s in ss]);
        Include(~ans,<v[1],v[2],v[3],v[4]>);
    end for;
    return Sort(Setseq(ans));
end function;

function ReadTSV(path)
    ans := {};
    lines := Split(Read(path),"\n");
    if #lines eq 0 then return ans; end if;
    head := Split(lines[1],"\t");
    pos := [Index(head,s) : s in ["primitive_a","primitive_b","primitive_c","primitive_d"]];
    if 0 in pos then return ans; end if;
    for line in lines[2..#lines] do
        if #line eq 0 then continue; end if;
        ss := Split(line,"\t");
        if #ss lt Max(pos) then continue; end if;
        v := Primitive([StringToInteger(ss[j]) : j in pos]);
        Include(~ans,<v[1],v[2],v[3],v[4]>);
    end for;
    return ans;
end function;

primes := [17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199];
bankset := Seqset(ReadBank("data/tor2228_bank.txt"));
for path in [
    "results/target_22224_offrectangle_new_curves_pairfiber_candidates.tsv",
    "results/target_22224_offrectangle_new_curves_pairfiber_small.tsv"
] do
    bankset join:= ReadTSV(path);
end for;
bank := Sort(Setseq(bankset));
tested := 0; offrect := 0; obstructed := 0; survivors := [];
print "TARGET_22224_A2228_TWIST_SIEVE_START","bank",#bank;
for tt in bank do
    v := [Z!tt[i] : i in [1..4]];
    if Rectangle(v) then continue; end if;
    offrect +:= 1;
    constraints := [];
    killed := false;
    for p in primes do
        k := GF(p); R<X> := PolynomialRing(k);
        fp := X*&*[X+(k!(z mod p))^2 : z in v];
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
        eta := k!2; while IsSquare(eta) do eta +:= 1; end while;
        np := #Jacobian(HyperellipticCurve(fp));
        nm := #Jacobian(HyperellipticCurve(eta*fp));
        plus := np mod 192 eq 0; minus := nm mod 192 eq 0;
        tested +:= 1;
        if not plus and not minus then
            obstructed +:= 1; killed := true;
            print "TWIST_BANK_OBSTRUCTION",v,"prime",p,
                  "plus_order",np,"minus_order",nm,
                  "v2",Valuation(np,2),Valuation(nm,2),
                  "v3",Valuation(np,3),Valuation(nm,3);
            break;
        elif plus ne minus then
            Append(~constraints,<p,plus select 1 else -1>);
        end if;
    end for;
    if not killed then
        Append(~survivors,<v,constraints>);
        print "TWIST_BANK_SURVIVOR",v,"constraints",constraints;
    end if;
end for;
print "TARGET_22224_A2228_TWIST_SIEVE_DONE","offrectangle",offrect,
      "local_tests",tested,"obstructed",obstructed,
      "survivors",#survivors;
quit;
