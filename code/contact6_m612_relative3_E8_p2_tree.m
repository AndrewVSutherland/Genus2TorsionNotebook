//////////////////////////////////////////////////////////////////////
// RIGOROUS 2-adic solution tree for E8's plane model
//     P(w,e) = e^4 * mp8(w) = e^4*w^8 + (216e^2+72e-24)e^2*w^4
//              + (-1296e^3-1728e^2-432e+64)e*w^2
//              + (-3888e^4-2592e^3+432e^2+288e-48)   in Z[w,e].
// Charts: (w,e), (1/W,e), (w,1/E), (1/W,1/E) with W,E in 2Z_2, cleared.
// Tree over Z_2^2: node = (w0,e0 mod 2^n); prune if v2(P(w0,e0)) < n;
// HENSEL ABORT if v2(P) > 2*min(v2(dP/dw), v2(dP/de)) (smooth lift =>
// a 2-adic point EXISTS).  All-trees-die => E8(Q_2) = empty on smooth
// affine charts; surviving deep nodes are reported for hand analysis.
// Usage: magma -b code/contact6_m612_relative3_E8_p2_tree.m
//////////////////////////////////////////////////////////////////////
SetColumnms := 0;
SetColumns(0); SetMemoryLimit(4*10^9);
Z := Integers();
R<w,ee> := PolynomialRing(Z, 2);
P := ee^4*w^8 + (216*ee^2+72*ee-24)*ee^2*w^4
     + (-1296*ee^3-1728*ee^2-432*ee+64)*ee*w^2
     + (-3888*ee^4-2592*ee^3+432*ee^2+288*ee-48);

// chart substitutions, denominators cleared minimally (unit scalings kept)
function ChartPolys(P, R)
    w := R.1; ee := R.2;
    charts := [];
    Append(~charts, <"w,e", P>);
    // w = 1/W (W in 2Z_2): multiply by W^8
    PW := ee^4 + (216*ee^2+72*ee-24)*ee^2*w^4 + (-1296*ee^3-1728*ee^2-432*ee+64)*ee*w^6
          + (-3888*ee^4-2592*ee^3+432*ee^2+288*ee-48)*w^8;
    Append(~charts, <"1/W,e", PW>);
    // e = 1/E (E in 2Z_2): multiply by E^4
    PE := w^8 + (216+72*ee-24*ee^2)*w^4 + (-1296-1728*ee-432*ee^2+64*ee^3)*w^2
          + (-3888-2592*ee+432*ee^2+288*ee^3-48*ee^4);
    Append(~charts, <"w,1/E", PE>);
    // both: clear W^8 E^4
    PWE := 1 + (216+72*ee-24*ee^2)*w^4 + (-1296-1728*ee-432*ee^2+64*ee^3)*w^6
           + (-3888-2592*ee+432*ee^2+288*ee^3-48*ee^4)*w^8;
    Append(~charts, <"1/W,1/E", PWE>);
    return charts;
end function;

// verify chart polys against symbolic substitution
FQ<a,b> := FunctionField(Rationals(), 2);
Pq := b^4*a^8 + (216*b^2+72*b-24)*b^2*a^4 + (-1296*b^3-1728*b^2-432*b+64)*b*a^2
      + (-3888*b^4-2592*b^3+432*b^2+288*b-48);
charts := ChartPolys(P, R);
chk2 := Evaluate(Pq, [1/a, b])*a^8;
chk3 := Evaluate(Pq, [a, 1/b])*b^4;
chk4 := Evaluate(Pq, [1/a, 1/b])*a^8*b^4;
sym := [Pq, chk2, chk3, chk4];
for i in [1..4] do
    lhs := &+[FQ | ]; // rebuild chart poly in FQ
    cp := charts[i][2];
    lhs := Evaluate(cp, [a, b]);
    assert lhs eq sym[i];
end for;
print "chart polynomials verified symbolically";

NMAX := 26;
procedure RunTree(name, Pc, R, restrictW, restrictE)
    // restrictW/E: true = variable must be in 2Z_2 (chart var = 1/orig)
    w := R.1; ee := R.2;
    Pw := Derivative(Pc, 1); Pe := Derivative(Pc, 2);
    nodes := [ <Z!0, Z!0, 0> ];   // (w0, e0, n): w=w0+2^n*Z2 etc.
    // apply restrictions at depth 1 later
    hensel := false; hw := <0,0,0>;
    for n in [0..NMAX-1] do
        newnodes := [];
        for nd in nodes do
            w0 := nd[1]; e0 := nd[2];
            for dw in [0,1] do for de in [0,1] do
                w1 := w0 + dw*2^n; e1 := e0 + de*2^n; n1 := n+1;
                if n eq 0 and restrictW and (w1 mod 2) eq 1 then continue; end if;
                if n eq 0 and restrictE and (e1 mod 2) eq 1 then continue; end if;
                val := Evaluate(Pc, [w1, e1]);
                vP := val eq 0 select 10^6 else Valuation(val, 2);
                if vP lt n1 then continue; end if;   // no lift in this branch
                // Hensel test
                vw := Evaluate(Pw, [w1, e1]); ve := Evaluate(Pe, [w1, e1]);
                vmin := Min(vw eq 0 select 10^6 else Valuation(vw,2),
                            ve eq 0 select 10^6 else Valuation(ve,2));
                if vP gt 2*vmin and vmin lt 10^6 then
                    hensel := true; hw := <w1, e1, n1>;
                    break nd;
                end if;
                Append(~newnodes, <w1, e1, n1>);
            end for; end for;
        end for;
        if hensel then break; end if;
        nodes := newnodes;
        if #nodes eq 0 then
            printf "chart %o : tree DIED at depth %o  => no smooth Q_2 point in chart\n", name, n+1;
            return;
        end if;
        if #nodes gt 4000 then
            printf "chart %o : node explosion (%o) at depth %o -- report and stop\n", name, #nodes, n+1;
            return;
        end if;
    end for;
    if hensel then
        printf "chart %o : HENSEL POINT at (w,e) = (%o, %o) mod 2^%o  => E8(Q_2) NONEMPTY\n",
            name, hw[1], hw[2], hw[3];
    else
        printf "chart %o : %o live nodes at depth %o (inconclusive; deepen or resolve)\n",
            name, #nodes, NMAX;
        for nd in nodes[1..Min(#nodes,10)] do
            printf "   live: w=%o e=%o mod 2^%o\n", nd[1], nd[2], nd[3];
        end for;
    end if;
end procedure;

for i in [1..4] do
    nm := charts[i][1];
    RunTree(nm, charts[i][2], R, i in [2,4], i in [3,4]);
end for;
print "E8_P2_TREE_DONE";
quit;
