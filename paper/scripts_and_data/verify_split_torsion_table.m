// verify_split_torsion_table.m -- verifies the exact torsion subgroup of
// every curve in Table 2 (tab:splitcensus) of the paper: for each row of
// table2.txt (the machine-readable table data from which table2.tex is
// typeset by make_tables.py), the displayed curve's Jacobian has rational
// torsion subgroup exactly the displayed group; when the row carries a
// factored display, the factorization is checked to multiply out to f.
// Run from this directory:
//   magma -b verify_split_torsion_table.m
// Torsion only: geometric splitness is certified separately by
// verify_split_certificates.m, which cross-checks its curves row-for-row
// against the same table2.txt.
SetColumns(0);
SetMemoryLimit(8*10^9);
R<x> := PolynomialRing(Rationals());

procedure CheckTorsion(I, fh)
    assert I eq Invariants(AbelianGroup(TorsionSubgroup(Jacobian(
        SimplifiedModel(HyperellipticCurve(R!fh[1], R!fh[2]))))));
end procedure;

t0 := Cputime();
n := 0;
for line in Split(Read("table2.txt"), "\n") do
    if #line eq 0 or line[1] eq "#" then continue; end if;
    parts := Split(line, "|");
    assert #parts eq 7;   // group | [[f],[h]] | label | display | source | route | comment
    I := eval parts[1];
    fh := eval parts[2];
    if parts[4] ne "-" then   // factored display: product must equal f, h = 0
        assert &*[R | R!fac : fac in eval parts[4]] eq R!fh[1] and fh[2] eq [0];
    end if;
    tag := parts[3] eq "-" select parts[7] else parts[3];
    printf "checking %-14o  (%o)\n", I, tag;
    CheckTorsion(I, fh);
    n +:= 1;
end for;
assert n eq 77;
printf "ALL %o CHECKS PASSED (%.1o s)\n", n, Cputime()-t0;
quit;
