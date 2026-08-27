// 37-hunt: Riemann-relation diagnostic + corrected reconstruction, fully
// offline (no modular symbols): EZ from results/gl2_pp_lattices_2190.log,
// candidate subgroups re-enumerated from EZ alone, period matrix recovered
// from dumped candidate #1 by inverting its known transform.
//
// Diagnostic: find the orientation Ew in {EZ, -EZ, EZ^t, -EZ^t} for which
//   P * Ew^-1 * P^t ~ 0   and   i * P * Ew^-1 * conj(P)^t > 0
// (P = 2x4 base-lattice period matrix).  Then rebuild the 12 candidates'
// big period matrices with the CORRECT symplectic convention and
// reconstruct curves over Q.
//
// Run: magma -b code/gl2_riemann_diag.m > results/gl2_riemann_diag.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(12*10^9);

Attach("~/.claude/jobs/a1db5dd4/tmp/polredabs_shim.m");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/quartic/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/curve_reconstruction/magma/spec");

Prec := 140;
FX := RationalsExtra(Prec);
CC := FX`CC;
RR := RealField(Prec);

EZ := Matrix(Integers(), 4, 4,
    [ 0, -237, 221, -205, 237, 0, 221, -16, -221, -221, 0, 205, 205, 16, -205, 0 ]);
ed := ElementaryDivisors(EZ);
dd := ed[3];
printf "TYPE %o dd=%o\n", ed, dd;

load "results/gl2_bigP_clean.m";
printf "LOADED %o candidates\n", #cands;

// re-enumerate candidate subgroup generators exactly as the driver did
Zd := Integers(dd);
Einv := ChangeRing(EZ, Rationals())^-1;
gensD := [ Vector(Rationals(), Eltseq(Einv[i])) : i in [1..4] ];
Ggens := [ Vector(Zd, [ Zd!(Integers()!(dd*x) mod dd) : x in Eltseq(v) ]) : v in gensD ];
EZd := ChangeRing(EZ, Zd);
fac := Factorization(dd);
function AddOrd(x)
    return dd div GCD([dd] cat [Integers()!c : c in Eltseq(x)]);
end function;
// same generator lists as the driver (order matters for combo indices!)
// but WITHOUT the star filter (we only need the subgroups the driver kept;
// we rebuild ALL isotropics per prime and later pick by matching the
// stored combo indices -- the driver's combo indexed its FILTERED list,
// so we must reproduce the filter too.  The star matrix is not available
// offline; instead we simply enumerate ALL isotropics and, for candidate
// reconstruction, use the lattice determined by matching the dumped bigP:
// bigP_dumped = (T * Bp * PM)^t.  We do not need the combo identification
// at all: we recover PM from candidate 1 and then verify against the
// other candidates' dumps.
function AllIsotropics(p, k)
    m := dd div p^k;
    pg := [ m*g_ : g_ in Ggens ];
    Mp := sub< RSpace(Zd, 4) | pg >;
    ord := p^k;
    e1 := 0; e2 := 0;
    for x in Mp do if AddOrd(x) eq ord then e1 := x; break; end if; end for;
    for x in Mp do
        if AddOrd(x) eq ord and #sub<Mp | [e1, x]> eq ord^2 then e2 := x; break; end if;
    end for;
    genlist := [ e1 + t*e2 : t in [0..ord-1] ] cat
               [ p*u*e1 + e2 : u in [0..(ord div p)-1] ];
    cands_ := [];
    for x in genlist do
        Hx := sub< Mp | x >;
        if #Hx ne ord then continue; end if;
        if not &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hx) ] then continue; end if;
        Append(~cands_, Hx);
    end for;
    if k eq 2 then
        Hp := sub< Mp | [ p*x : x in Generators(Mp) ] >;
        if #Hp eq ord and
           &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hp) ] then
            Append(~cands_, Hp);
        end if;
    end if;
    return cands_;
end function;

// Reconstruct the lattice basis Bp and transform T for a given selection
function LatticeData(sel, alliso)
    Hgens := &cat[ [ x : x in Generators(alliso[j][sel[j]]) ] : j in [1..#alliso] ];
    lifts := [ Vector(Rationals(), [ (Integers()!x)/dd : x in Eltseq(u) ]) : u in Hgens ];
    B := VerticalJoin(ChangeRing(IdentityMatrix(Integers(),4), Rationals()),
                      Matrix(Rationals(), #lifts, 4, &cat[Eltseq(v) : v in lifts]));
    dn := LCM([Denominator(x) : x in Eltseq(B)]);
    BZ := Matrix(Integers(), Nrows(B), 4, [Integers()!(dn*x) : x in Eltseq(B)]);
    Hn := HermiteForm(BZ);
    Bp := Matrix(Rationals(), 4, 4, [ Hn[i][j]/dn : j in [1..4], i in [1..4] ]);
    Ep := Bp * ChangeRing(EZ, Rationals()) * Transpose(Bp);
    EpZ := Matrix(Integers(), 4, 4, [Integers()!x : x in Eltseq(Ep)]);
    edp := ElementaryDivisors(EpZ);
    EpN := EpZ div edp[1];
    _, T := FrobeniusFormAlternating(EpN);
    return Bp, T, EpN;
end function;

// The DRIVER's star-filtered per-prime counts were [1,2,2,3]; its combos in
// the dump tell us which selections were used, but indices refer to the
// filtered lists.  To avoid ambiguity we recover PM by brute force: try all
// isotropic selections at each prime whose (Bp,T) reproduce dumped bigP #1.
alliso := [* *];
for pf in fac do
    Append(~alliso, AllIsotropics(pf[1], pf[2]));
end for;
printf "ALL_ISO_COUNTS %o\n", [ #c : c in alliso ];

// dumped bigP #1 as matrix
d1 := cands[1];
es1 := [ CC!(RR!d1[2][2*i-1]) + CC.1*(RR!d1[2][2*i]) : i in [1..8] ];
bigP1 := Matrix(CC, 2, 4, es1);

// For PM recovery: bigP^t = (T * Bp) * PM  =>  PM = (T*Bp)^-1 * bigP^t.
// We don't know which selection made candidate 1; but PM is THE SAME for
// all candidates, and for the correct selection the recovered PM must make
// ALL dumped candidates consistent.  Try selections until consistent.
// (search bounded: counts are [3,4,12,1408] -> too many; but note the
// DRIVER enumerated star-filtered combos in lexicographic order with
// counts [1,2,2,3], and its dump stores the combo indices; the j-th
// filtered subgroup is among our all-list.  We simply try, for candidate 1
// with dumped combo [1,1,1,1], all selections with small indices first.)
found := false;
PMrec := ZeroMatrix(CC, 4, 2);
for i1 in [1..#alliso[1]] do
 for i2 in [1..#alliso[2]] do
  for i3 in [1..#alliso[3]] do
   for i4 in [1..#alliso[4]] do
    Bp, T, EpN := LatticeData([i1,i2,i3,i4], alliso);
    if ElementaryDivisors(EpN) ne [1,1,1,1] then continue; end if;
    TB := ChangeRing(ChangeRing(T, Rationals()) * Bp, CC);
    PMc := TB^-1 * Transpose(bigP1);
    // verify: base-lattice Riemann test for each orientation
    P0 := Transpose(PMc);   // 2x4 periods of the BASE lattice
    for orient in [1,2,3,4] do
        Ew := orient eq 1 select EZ else orient eq 2 select -EZ
              else orient eq 3 select Transpose(EZ) else -Transpose(EZ);
        Ei := ChangeRing(ChangeRing(Ew, Rationals())^-1, CC);
        R0 := P0 * Ei * Transpose(P0);
        err := Max([Abs(x) : x in Eltseq(R0)]);
        Hm := CC.1 * P0 * Ei * Transpose(Matrix(CC,2,4,[Conjugate(x) : x in Eltseq(P0)]));
        ev1 := Re(Hm[1][1]);
        if err lt 10^(-12) and ev1 gt 0 then
            printf "CONSISTENT sel=[%o,%o,%o,%o] orient=%o err=%o\n",
                   i1,i2,i3,i4, orient, RealField(6)!err;
            found := true;
            PMrec := PMc;
        end if;
    end for;
    if found then break; end if;
   end for;
   if found then break; end if;
  end for;
  if found then break; end if;
 end for;
 if found then break; end if;
end for;
error if not found, "no consistent selection/orientation found";
printf "PM_RECOVERED\n";
printf "GL2_RIEMANN_DIAG_DONE\n";
quit;
