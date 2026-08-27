// Direct computation of the cuspidal image in the optimal quotient A_f,
// via the rational period mapping on modular symbols (Manin-Drinfeld).
//
// For f = 2190.2.a.v (target 37) / 1830.2.a.q (target 31):
//   phi := RationalMapping(Af)  maps  M_2(Gamma_0(N);Q) -> Q^4,
//   Lambda := phi(H_1(X_0(N),Z))  is the period lattice of the OPTIMAL
//   QUOTIENT torus; for each cusp c (all rational: N squarefree), the
//   image of the class [(c)-(oo)] in A_f(Q) is phi({c,oo}) mod Lambda,
//   rational by Galois-equivariance and torsion by Manin-Drinfeld.
//   The subgroup order = index [ Lambda + sum_c Z*phi({c,oo}) : Lambda ].
// 37 | index  =>  THEOREM: rational 37-torsion point on A_f.
//
// Run: magma -b Lv:=2190 code/gl2_cusp_image_order.m > results/gl2_cuspimg_2190.log

SetColumns(0);
SetSeed(1);
if not assigned Lv then Lv := 2190; elif Type(Lv) eq MonStgElt then Lv := StringToInteger(Lv); end if;
if not assigned MemGB then MemGB := 16; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

if Lv eq 2190 then
    trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
    targetprime := 37;
elif Lv eq 1830 then
    trtargets := [<7,0>, <11,8>, <13,-2>, <17,6>, <19,-6>];
    targetprime := 31;
elif Lv eq 23 then
    trtargets := [];   // single newform piece; validation: cusp 0 must have order 11
    targetprime := 11;
else
    error "unknown level";
end if;

printf "LEVEL %o TARGET %o\n", Lv, targetprime;
M := ModularSymbols(Lv, 2, 0);
S := CuspidalSubspace(M);
NS := NewSubspace(S);
D := NewformDecomposition(NS);
target := 0;
for i in [1..#D] do
    if Dimension(D[i]) ne 4 then continue; end if;
    ok := true;
    for tt in trtargets do
        if Trace(HeckeOperator(D[i], tt[1])) ne 2*tt[2] then ok := false; break; end if;
    end for;
    if ok then target := i; break; end if;
end for;
error if target eq 0, "piece not found";
Af := D[target];
printf "PIECE %o dim=%o\n", target, Dimension(Af);

phi := RationalMapping(Af);
// period lattice of the optimal quotient: image of integral cuspidal homology
IB := IntegralBasis(S);
imgs := [ phi(M!b) : b in IB ];
V := Universe(imgs);
n := Degree(V);
printf "PERIOD_SPACE_DIM %o LATTICE_GENS %o\n", n, #imgs;
LM0 := Matrix(Rationals(), [ Eltseq(v) : v in imgs ]);
den := LCM([ Denominator(x) : x in Eltseq(LM0) ]);
Zmat := Matrix(Integers(), [[ Integers()!(den*x) : x in Eltseq(v) ] : v in imgs]);
H := HermiteForm(Zmat);
Hrows := [ H[i] : i in [1..Nrows(H)] | not IsZero(H[i]) ];
error if #Hrows ne n, "Hermite rank mismatch";
LamZ := Matrix(Integers(), #Hrows, n, &cat[ Eltseq(r) : r in Hrows ]);   // basis of den*Lambda

// cusp images
cusps := [ c : c in Cusps(Gamma0(Lv)) ];
printf "N_CUSPS %o\n", #cusps;
orders := [];
for c in cusps do
    cc := Eltseq(c);
    alpha := cc[2] eq 0 select Infinity() else cc[1]/cc[2];
    if alpha cmpeq Infinity() then continue; end if;
    x := M ! <1, [alpha, Infinity()]>;
    v := phi(x);
    vz := Vector(Rationals(), [ den*t : t in Eltseq(v) ]);
    // order of v mod Lambda: smallest k with k*vz in Z-span(LamZ)
    sol := Solution(ChangeRing(LamZ, Rationals()), vz);
    k := LCM([ Denominator(t) : t in Eltseq(sol) ]);
    Append(~orders, k);
    printf "CUSP %o image_order=%o\n", alpha, k;
end for;
printf "CUSP_IMAGE_ORDERS %o\n", orders;
ord_lcm := LCM(orders);
printf "SUBGROUP_ORDER_LCM_LOWER %o\n", ord_lcm;
if ord_lcm mod targetprime eq 0 then
    printf "THEOREM: A_f (optimal quotient) has a rational point of order %o (cuspidal image)\n", targetprime;
else
    printf "CUSP_IMAGES_INSUFFICIENT lcm=%o\n", ord_lcm;
end if;
printf "GL2_CUSPIMG_DONE level=%o\n", Lv;
quit;
