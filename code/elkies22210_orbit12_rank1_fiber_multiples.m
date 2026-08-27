/*
Enumerate Mordell--Weil multiples on rational fibers of the S_3 quotient
of the first orbit-12 radicand cover.

On a fixed value of

  s = r3*r4+r3*r5+r4*r5,

the quotient curve is

  v^2 = q*(s-q)*(2-q^2+(q+2)*s)

and has Weierstrass model

  Y^2 = X^3+(s^2-2*s-2)*X^2-4*s^2*(s+1)*X
          +4*s^2*(s+1)^2,
  q = 2*s*(s+1)/X.

The quotient forgets the ordering of the three complementary roots.  A
candidate lifts to the labelled Clebsch--Klein surface exactly when

  z^3+(1+q)z^2+s*z-(1+q)*(q-s)

splits completely over Q.  Every such lift is checked against the two CK
identities, smoothness, and all four literal Stoll--Zarhin radicands.

Usage:

  magma -b integer_bound:=30 s_height:=0 multiple_bound:=100 \
    code/elkies22210_orbit12_rank1_fiber_multiples.m

Here integer_bound searches all integral s in the indicated symmetric
interval.  If s_height is positive, reduced nonintegral s=a/b with
|a|,b <= s_height are added.  The known split-cubic regression fiber
s=59/49 is always included unless include_seed:=0.
*/

if not assigned multiple_bound then
    multiple_bound := 100;
elif Type(multiple_bound) eq MonStgElt then
    multiple_bound := StringToInteger(multiple_bound);
end if;
require multiple_bound ge 1 : "multiple_bound must be positive";

if not assigned integer_bound then
    integer_bound := 30;
elif Type(integer_bound) eq MonStgElt then
    integer_bound := StringToInteger(integer_bound);
end if;
if not assigned s_height then
    s_height := 0;
elif Type(s_height) eq MonStgElt then
    s_height := StringToInteger(s_height);
end if;
if not assigned include_seed then
    include_seed := 1;
elif Type(include_seed) eq MonStgElt then
    include_seed := StringToInteger(include_seed);
end if;
require integer_bound ge 0 and s_height ge 0 : "bounds must be nonnegative";

Q := Rationals();
P<z> := PolynomialRing(Q);
SetSeed(1);
cpu_start := Cputime();
s_values := [Q!s : s in [-integer_bound..integer_bound] |
             s notin [0,-1]];
if s_height gt 0 then
    for b in [2..s_height] do
        for a in [-s_height..s_height] do
            if a ne 0 and a ne -b and GCD(Abs(a),b) eq 1 then
                Append(~s_values,Q!(a/b));
            end if;
        end for;
    end for;
end if;
if include_seed ne 0 then Append(~s_values,Q!(59/49)); end if;
s_values := Setseq(Seqset(s_values));
Sort(~s_values);

function LiteralRadicands(rr)
    aa := [x^2 : x in rr];
    return [
        -(aa[1]-aa[3])*(aa[1]-aa[4])*(aa[1]-aa[5]),
         (aa[3]-aa[2])*(aa[1]-aa[4])*(aa[1]-aa[5]),
         (aa[4]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[5]),
         (aa[5]-aa[2])*(aa[1]-aa[3])*(aa[1]-aa[4])
    ];
end function;

function SmoothCK(rr)
    if &or[x eq 0 : x in rr] then return false; end if;
    if &+[x : x in rr] ne 0 or &+[x^3 : x in rr] ne 0 then
        return false;
    end if;
    aa := [x^2 : x in rr];
    return &and[aa[i] ne aa[j] : i,j in [1..5] | i lt j];
end function;

total_points := 0;
finite_q := 0;
split_one_root := 0;
split_three_roots := 0;
smooth_lifts := 0;
first_square_failures := 0;
full_hits := 0;
max_mask_weight := 0;
rank_zero_fibers := 0;
positive_rank_fibers := 0;
ambiguous_rank_fibers := 0;

printf "ELKIES22210_ORBIT12_RANK1_FIBER_MULTIPLES\n";
printf "multiple_bound %o\n", multiple_bound;
printf "integer_bound %o\n", integer_bound;
printf "s_height %o\n", s_height;
printf "fiber_count %o\n", #s_values;

/* Regression: this point is a genuine smooth split lift of the first
   radicand quotient, but only G0 is a square. */
sseed := Q!(59/49);
qseed := Q!(8/7);
vseed := Q!(192/343);
Xseed := 2*sseed*(sseed+1)/qseed;
Yseed := 2*sseed*(sseed+1)*vseed/qseed^2;
Eseed := EllipticCurve([Q | 0, sseed^2-2*sseed-2, 0,
                       -4*sseed^2*(sseed+1),
                        4*sseed^2*(sseed+1)^2]);
eseed := Eseed![Xseed,Yseed,1];
assert qseed eq 2*sseed*(sseed+1)/eseed[1];
cseed := z^3+(1+qseed)*z^2+sseed*z-(1+qseed)*(qseed-sseed);
rseed := [Q!-9/7,Q!-5/7,Q!-1/7];
assert &and[Evaluate(cseed,a) eq 0 : a in rseed];
Gseed := LiteralRadicands([Q!1,qseed] cat rseed);
seedmask := [IsSquare(g) : g in Gseed];
assert seedmask eq [true,false,false,false];
printf "REGRESSION s=%o q=%o Epoint=%o roots=%o mask=%o\n",
       sseed,qseed,eseed,rseed,seedmask;

for s in s_values do
    E := EllipticCurve([Q | 0, s^2-2*s-2, 0,
                       -4*s^2*(s+1), 4*s^2*(s+1)^2]);
    lo, hi := RankBounds(E);
    if lo eq 0 then
        if hi eq 0 then rank_zero_fibers +:= 1;
        else ambiguous_rank_fibers +:= 1;
        end if;
        printf "FIBER_SKIP s=%o rank_bounds=%o,%o\n",s,lo,hi;
        continue;
    end if;
    positive_rank_fibers +:= 1;

    gens := Generators(E);
    free := [g : g in gens | Order(g) eq 0];
    assert #free ge 1;

    T, tmap := TorsionSubgroup(E);
    torsion_points := [tmap(a) : a in T];
    seen_q := {};

    fiber_points := 0;
    fiber_finite_q := 0;
    fiber_one_root := 0;
    fiber_split := 0;
    fiber_smooth := 0;
    fiber_hits := 0;

    for gi in [1..#free] do
      gen := free[gi];
      for n in [1..multiple_bound] do
            for tt in torsion_points do
                ep := n*gen + tt;
                fiber_points +:= 1;
                total_points +:= 1;
                if IsIdentity(ep) or ep[1] eq 0 then continue; end if;

                q := Q!(2*s*(s+1)/ep[1]);
                if q in seen_q then continue; end if;
                Include(~seen_q,q);
                fiber_finite_q +:= 1;
                finite_q +:= 1;

                cubic := z^3+(1+q)*z^2+s*z-(1+q)*(q-s);
                fac := Factorization(cubic);
                linear := [row[1] : row in fac | Degree(row[1]) eq 1];
                if #linear gt 0 then
                    fiber_one_root +:= 1;
                    split_one_root +:= 1;
                    printf "REDUCIBLE_CUBIC s=%o gen=%o n=%o q=%o factorization=%o\n",
                           s,gi,n,q,fac;
                end if;
                if #fac ne 3 or not &and[Degree(row[1]) eq 1 and
                                         row[2] eq 1 : row in fac] then
                    continue;
                end if;

                roots := [ -Coefficient(row[1],0)/Coefficient(row[1],1)
                           : row in fac ];
                assert #Seqset(roots) eq 3;
                fiber_split +:= 1;
                split_three_roots +:= 1;

                rr := [Q!1,q] cat roots;
                if not SmoothCK(rr) then continue; end if;
                fiber_smooth +:= 1;
                smooth_lifts +:= 1;

                GG := LiteralRadicands(rr);
                sq := [IsSquare(g) : g in GG];
                if not sq[1] then first_square_failures +:= 1; end if;
                weight := #[b : b in sq | b];
                max_mask_weight := Max(max_mask_weight,weight);
                printf "SPLIT_LIFT s=%o gen=%o n=%o q=%o roots=%o mask=%o G=%o\n",
                       s,gi,n,q,roots,sq,GG;
                if &and sq then
                    fiber_hits +:= 1;
                    full_hits +:= 1;
                    printf "HIT r=%o\n", rr;
                end if;
            end for;
      end for;
    end for;

    printf "FIBER s=%o rank_bounds=%o,%o torsion=%o free_generators=%o points=%o finite_q=%o one_rational_root=%o split=%o smooth=%o hits=%o\n",
           s,lo,hi,Invariants(T),free,fiber_points,fiber_finite_q,
           fiber_one_root,fiber_split,fiber_smooth,fiber_hits;
end for;

printf "rank_zero_fibers %o\n", rank_zero_fibers;
printf "ambiguous_rank_fibers %o\n", ambiguous_rank_fibers;
printf "positive_rank_fibers %o\n", positive_rank_fibers;
printf "total_points %o\n", total_points;
printf "finite_q %o\n", finite_q;
printf "cubics_with_rational_root %o\n", split_one_root;
printf "cubics_split_completely %o\n", split_three_roots;
printf "smooth_split_lifts %o\n", smooth_lifts;
printf "first_square_failures %o\n", first_square_failures;
printf "max_exact_square_mask_weight %o\n", max_mask_weight;
printf "full_hits %o\n", full_hits;
printf "cpu_seconds %o\n", Cputime(cpu_start);
printf "DONE\n";

quit;
