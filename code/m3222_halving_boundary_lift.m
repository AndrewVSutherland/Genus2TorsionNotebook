//////////////////////////////////////////////////////////////////////
//  Hensel lifting on the boundary charts for the symmetric halving
//  cover in M_1(8,2,2).
//
//  Variables:
//      s,p,A,B,C,D,E,z
//
//  Equations:
//      f_{s,p}(X) - L*(X+1)*(X^2+A*X+B)^2 = (C*X^2+D*X+E)^2,
//      z^2 = s^2 - 4p.
//
//  We start from split boundary solutions modulo ell and lift by solving
//  the linearized equations from ell^k to ell^(k+1).
//
//  Typical runs from torsion_jac:
//      magma -b ell:=7 levels:=3 code/m3222_halving_boundary_lift.m
//      magma -b ell:=11 levels:=3 code/m3222_halving_boundary_lift.m
//////////////////////////////////////////////////////////////////////

if not assigned ell then
    ell := 7;
elif Type(ell) eq MonStgElt then
    ell := StringToInteger(ell);
end if;
if not assigned levels then
    levels := 3;
elif Type(levels) eq MonStgElt then
    levels := StringToInteger(levels);
end if;
if not assigned max_solutions then
    max_solutions := 500000;
elif Type(max_solutions) eq MonStgElt then
    max_solutions := StringToInteger(max_solutions);
end if;
if not assigned boundary_only then
    boundary_only := true;
elif Type(boundary_only) eq MonStgElt then
    boundary_only := boundary_only in {"true", "True", "1", "yes"};
end if;

Z := Integers();
R<s,pp,A,B,C,D,E,z> := PolynomialRing(Z, 8);
PX<X> := PolynomialRing(R);

qtilde := -X^2 + (pp*s - s^2 + 2*pp - s - 2)*X - (s^2 - pp + s + 1);
f := ((pp-s+1)*X^2 + (2-s)*X + 1)*((s+2)*X + 1)*qtilde;
L := Coefficient(f, 5);
a := X^2 + A*X + B;
rootpoly := C*X^2 + D*X + E;
res := f - L*(X+1)*a^2 - rootpoly^2;
eqs := [Coefficient(res, i) : i in [0..4]];
Append(~eqs, z^2 - (s^2 - 4*pp));
jac := [[Derivative(eqs[i], j) : j in [1..8]] : i in [1..#eqs]];

labels := [
    "p=0",
    "s+1=0",
    "s+2=0",
    "p-s+1=0",
    "s-p+1=0",
    "delta=0",
    "qdisc=0",
    "collision=0"
];

function BoundaryValues(ss, p0)
    return [
        p0,
        ss + 1,
        ss + 2,
        p0 - ss + 1,
        ss - p0 + 1,
        ss^2 - 4*p0,
        2*ss^2 + 3*ss + p0 + 1,
        ss^3 - ss^2*p0 + ss^2 - 4*ss*p0 - 4*p0
    ];
end function;

function LabelString(ss, p0)
    vals := BoundaryValues(ss, p0);
    out := "";
    for i in [1..#labels] do
        if vals[i] eq 0 then
            if #out gt 0 then
                out cat:= ",";
            end if;
            out cat:= labels[i];
        end if;
    end for;
    if #out eq 0 then
        return "open";
    end if;
    return out;
end function;

function EvalInt(poly, vals)
    return Z!Evaluate(poly, vals);
end function;

function Key(vals, modulus)
    return &cat [Sprint(vals[i] mod modulus) cat ":" : i in [1..#vals]];
end function;

procedure AppendUnique(~out, ~seen, entry, modulus)
    key := Key(entry[1], modulus);
    if key notin seen then
        Include(~seen, key);
        Append(~out, entry);
    end if;
end procedure;

function InitialSolutions(ell)
    F := GF(ell);
    PF<T> := PolynomialRing(F);
    PFX<XF> := PolynomialRing(F);
    out := [];
    seen := {};

    for ss in F do
        for p0 in F do
            bvals := BoundaryValues(ss, p0);
            is_boundary := &or [v eq 0 : v in bvals];
            if boundary_only and not is_boundary then
                continue;
            end if;

            zroots := Roots(T^2 - (ss^2 - 4*p0));
            if #zroots eq 0 then
                continue;
            end if;

            qtildeF := -XF^2 + (p0*ss - ss^2 + 2*p0 - ss - 2)*XF
                       - (ss^2 - p0 + ss + 1);
            fF := ((p0-ss+1)*XF^2 + (2-ss)*XF + 1)*((ss+2)*XF + 1)*qtildeF;
            LF := Coefficient(fF, 5);

            for AA in F do
                for BB in F do
                    aa := XF^2 + AA*XF + BB;
                    residual := fF - LF*(XF+1)*aa^2;
                    square_ok, rr := IsSquare(residual);
                    if not square_ok then
                        continue;
                    end if;

                    roots := [rr];
                    if rr ne 0 then
                        Append(~roots, -rr);
                    end if;

                    for rpoly in roots do
                        CC := Coefficient(rpoly, 2);
                        DD := Coefficient(rpoly, 1);
                        EE := Coefficient(rpoly, 0);
                        for zr in zroots do
                            vals := [
                                Z!ss, Z!p0, Z!AA, Z!BB,
                                Z!CC, Z!DD, Z!EE, Z!zr[1]
                            ];
                            label := LabelString(ss, p0);
                            basekey := (Z!ss) + ell*(Z!p0);
                            AppendUnique(~out, ~seen, <vals, basekey, label>, ell);
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;

    return out;
end function;

function LiftEntry(entry, modulus, ell)
    vals := entry[1];
    basekey := entry[2];
    label := entry[3];
    F := GF(ell);

    rows := [];
    rhs := [];
    for i in [1..#eqs] do
        val := EvalInt(eqs[i], vals);
        if val mod modulus ne 0 then
            print "WARNING: entry not a solution modulo", modulus, "eq", i, "val", val;
            return [];
        end if;
        Append(~rhs, F!((-(val div modulus)) mod ell));

        row := [];
        for j in [1..8] do
            Append(~row, F!(EvalInt(jac[i][j], vals) mod ell));
        end for;
        Append(~rows, row);
    end for;

    M := Matrix(F, rows);
    b := Vector(F, rhs);
    MSolve := Transpose(M);
    if not IsConsistent(MSolve, b) then
        return [];
    end if;

    part := Solution(MSolve, b);
    K := Nullspace(MSolve);
    bas := Basis(K);
    dim := #bas;

    lifts := [];
    if dim eq 0 then
        tvecs := [part];
    else
        tvecs := [];
        for coeffs in CartesianPower(F, dim) do
            tv := part;
            for i in [1..dim] do
                tv +:= coeffs[i]*bas[i];
            end for;
            Append(~tvecs, tv);
        end for;
    end if;

    for tv in tvecs do
        newvals := [vals[i] + modulus*(Z!tv[i]) : i in [1..8]];
        Append(~lifts, <newvals, basekey, label, dim>);
    end for;
    return lifts;
end function;

function BoundaryIntValuesFromEntry(entry)
    vals := entry[1];
    ss := vals[1];
    p0 := vals[2];
    return [
        p0,
        ss + 1,
        ss + 2,
        p0 - ss + 1,
        ss - p0 + 1,
        ss^2 - 4*p0,
        2*ss^2 + 3*ss + p0 + 1,
        ss^3 - ss^2*p0 + ss^2 - 4*ss*p0 - 4*p0
    ];
end function;

function ComboString(idxs)
    if #idxs eq 0 then
        return "none";
    end if;
    out := "";
    for n in [1..#idxs] do
        if n gt 1 then
            out cat:= ",";
        end if;
        out cat:= labels[idxs[n]];
    end for;
    return out;
end function;

procedure SummarizeBoundaryDepth(level, modulus, entries, ell)
    if level lt 2 then
        return;
    end if;
    stuck_counts := AssociativeArray();
    base_counts := AssociativeArray(Integers());
    base_none_counts := AssociativeArray(Integers());
    base_labels := AssociativeArray(Integers());
    generic_depth1 := 0;
    for entry in entries do
        bvals := BoundaryIntValuesFromEntry(entry);
        active := [i : i in [1..#labels] | bvals[i] mod ell eq 0];
        stuck := [i : i in [1..#labels] | bvals[i] mod modulus eq 0];
        combo := ComboString(stuck);
        bk := entry[2];
        if not IsDefined(base_counts, bk) then
            base_counts[bk] := 0;
            base_none_counts[bk] := 0;
            base_labels[bk] := entry[3];
        end if;
        base_counts[bk] +:= 1;
        if combo eq "none" then
            base_none_counts[bk] +:= 1;
        end if;
        if not IsDefined(stuck_counts, combo) then
            stuck_counts[combo] := 0;
        end if;
        stuck_counts[combo] +:= 1;
        if #active eq 1 and #stuck eq 0 then
            generic_depth1 +:= 1;
        end if;
    end for;
    print "  boundary_depth_summary_mod", modulus, "generic_single_depth1", generic_depth1;
    for combo in Keys(stuck_counts) do
        print "    still_zero_mod_level", combo, "count", stuck_counts[combo];
    end for;
    print "  base_residue_depth_summary";
    for bk in Keys(base_counts) do
        ss0 := bk mod ell;
        p0 := (bk div ell) mod ell;
        print "    base", <ss0,p0>, "label", base_labels[bk],
              "total", base_counts[bk], "none", base_none_counts[bk];
    end for;
end procedure;

procedure Summarize(level, modulus, entries)
    base_set := { Z | };
    label_counts := AssociativeArray();
    for entry in entries do
        Include(~base_set, entry[2]);
        lab := entry[3];
        if not IsDefined(label_counts, lab) then
            label_counts[lab] := 0;
        end if;
        label_counts[lab] +:= 1;
    end for;

    print "level", level, "modulus", modulus,
          "solutions", #entries, "base_residues_mod_ell", #base_set;
    for lab in Keys(label_counts) do
        print "  label", lab, "count", label_counts[lab];
    end for;
    SummarizeBoundaryDepth(level, modulus, entries, ell);
end procedure;

print "Boundary Hensel lifting for M_1(8,2,2) halving cover";
print "ell", ell, "levels", levels, "boundary_only", boundary_only,
      "max_solutions", max_solutions;

entries := InitialSolutions(ell);
modulus := ell;
Summarize(1, modulus, entries);

for lev in [2..levels] do
    new_entries := [];
    seen := {};
    dim_counts := AssociativeArray(Integers());
    dead := 0;
    for entry in entries do
        lifts := LiftEntry(entry, modulus, ell);
        if #lifts eq 0 then
            dead +:= 1;
            continue;
        end if;
        dim := lifts[1][4];
        if not IsDefined(dim_counts, dim) then
            dim_counts[dim] := 0;
        end if;
        dim_counts[dim] +:= 1;

        for lift in lifts do
            AppendUnique(~new_entries, ~seen, <lift[1], lift[2], lift[3]>, modulus*ell);
            if #new_entries gt max_solutions then
                print "Exceeded max_solutions while lifting to level", lev;
                print "partial solutions", #new_entries;
                quit;
            end if;
        end for;
    end for;

    print "lift_to_level", lev, "dead_previous_solutions", dead;
    for d in Keys(dim_counts) do
        print "  tangent_dim", d, "previous_solution_count", dim_counts[d];
    end for;

    entries := new_entries;
    modulus *:= ell;
    Summarize(lev, modulus, entries);
end for;

print "sample lifted solutions";
for i in [1..Minimum(#entries, 20)] do
    print entries[i];
end for;

quit;
