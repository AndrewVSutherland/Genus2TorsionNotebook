//////////////////////////////////////////////////////////////////////
// Finite-group scout for transverse lines through the exact split Z/60
// seed.  For G_a= a+x and F_t=F_0+t*G_a, a rational point of order 60
// injects into J(F_p) at every good p not dividing 60.  Consequently the
// exponent of J_t(F_p) must be divisible by 60.  Singular fibers are kept
// as boundary residues and are never used to exclude a rational t.
//
// This is a direction-ranking diagnostic.  The selected direction is
// recomputed by z60_split_seed_finite_masks.m before any height sieve.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

Z := Integers();
if not assigned amin then amin := -4; end if;
if not assigned amax then amax := 4; end if;
if Type(amin) eq MonStgElt then amin := StringToInteger(amin); end if;
if Type(amax) eq MonStgElt then amax := StringToInteger(amax); end if;

primes := [7,11,13,17,19,23,29,31];

function HasElementOfOrder60(inv)
    return exists{n : n in inv | (Z!n) mod 60 eq 0};
end function;

print "Z60_SPLIT_SEED_DIRECTION_SCOUT";
print "family F_t=F_0+t*(a+x), a in", [amin..amax];

for aa in [amin..amax] do
    score_num := 1;
    score_den := 1;
    records := [* *];
    for p in primes do
        k := GF(p); P<x> := PolynomialRing(k);
        f0 := -46250000*x^6 + 1761500625*x^4
              - 22332312000*x^2 + 94277468160;
        allowed := []; bad := [];
        for tt in k do
            f := f0 + tt*(k!aa+x);
            if Degree(f) notin {5,6} or Discriminant(f) eq 0 then
                Append(~bad,Z!tt);
                continue;
            end if;
            try
                A, mp := AbelianGroup(Jacobian(HyperellipticCurve(f)));
                inv := Invariants(A);
                if HasElementOfOrder60(inv) then
                    Append(~allowed,Z!tt);
                end if;
            catch e
                Append(~bad,Z!tt);
            end try;
        end for;
        // Include one projective infinity residue in the density score.
        score_num *:= #Seqset(allowed cat bad)+1;
        score_den *:= p+1;
        Append(~records,<p,Sort(Setseq(Seqset(allowed))),
                          Sort(Setseq(Seqset(bad)))>);
    end for;
    print "a",aa,"score",RealField(8)!score_num/score_den,
          "records",records;
end for;

quit;
