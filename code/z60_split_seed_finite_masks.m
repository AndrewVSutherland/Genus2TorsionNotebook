//////////////////////////////////////////////////////////////////////
// Rigorous necessary finite-group masks on the transverse line
//
//   F_t = -46250000*x^6 + 1761500625*x^4
//         -22332312000*x^2 + 94277468160 + t*x.
//
// At a prime p not dividing 60, a rational point of order 60 injects
// into the good reduction.  Hence the exponent of J_t(F_p) is divisible
// by 60.  Singular fibers and computation failures are retained in the
// bad set, so the output is a rigorous necessary mask.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

if not assigned prime_bound then prime_bound := 181; end if;
if Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned output_file then
    output_file := "data/z60_split_seed_tx_masks_p181.txt";
end if;

Z := Integers();
out := Open(output_file,"w");

function HasElementOfOrder60(inv)
    return exists{n : n in inv | (Z!n) mod 60 eq 0};
end function;

fprintf out,"# Necessary masks for F_t=F_0+t*x; singular residues retained.\n";
fprintf out,"# p : allowed ; bad ; allowed invariant-factor records\n";
print "Z60_SPLIT_SEED_TX_FINITE_MASKS", "prime_bound", prime_bound;

singleton_primes := [];
singleton_product := Z!1;
for p in [q : q in [7..prime_bound] |
               IsPrime(q) and q notin {2,3,5,37}] do
    k := GF(p); P<x> := PolynomialRing(k);
    f0 := -46250000*x^6 + 1761500625*x^4
          - 22332312000*x^2 + 94277468160;
    allowed := []; bad := []; records := [* *];
    for tt in k do
        f := f0 + tt*x;
        if Degree(f) notin {5,6} or Discriminant(f) eq 0 then
            Append(~bad,Z!tt);
            continue;
        end if;
        try
            A, mp := AbelianGroup(Jacobian(HyperellipticCurve(f)));
            inv := Invariants(A);
            if HasElementOfOrder60(inv) then
                Append(~allowed,Z!tt);
                Append(~records,<Z!tt,inv>);
            end if;
        catch e
            Append(~bad,Z!tt);
        end try;
    end for;
    allowed := Sort(Setseq(Seqset(allowed)));
    bad := Sort(Setseq(Seqset(bad)));
    if allowed eq [0] and #bad eq 0 then
        Append(~singleton_primes,p);
        singleton_product *:= p;
    end if;
    fprintf out,"%o : %o ; %o ; %o\n",p,allowed,bad,records;
    print "p",p,"allowed",allowed,"bad",bad;
end for;

fprintf out,"# singleton_primes %o\n",singleton_primes;
fprintf out,"# singleton_product %o\n",singleton_product;
print "singleton_primes",singleton_primes;
print "singleton_product",singleton_product;
print "wrote",output_file;

delete out;
quit;
