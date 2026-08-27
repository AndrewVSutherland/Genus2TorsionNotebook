//////////////////////////////////////////////////////////////////////
//  Dump finite open residue classes where the distinguished order-8
//  class in M_1(8,2,2) is divisible by 2.
//
//  The output is a simple text format:
//      p key key key ...
//  where key = u_res + p*v_res.
//////////////////////////////////////////////////////////////////////

if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned output_file then
    output_file := "data/m3222_halving_allowed_residues_p43.txt";
end if;

Z := Integers();
prime_list := [p : p in [7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73] |
                    p le prime_bound];

function IsDivisibleBy2Finite(J, D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2, invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function FiniteOpenData(u, v, F, PF)
    XF := PF.1;
    if u eq 0 or v eq 0 or u eq v then
        return false, _, _, _;
    end if;
    if u eq 1 or v eq 1 or u+v+1 eq 0 or u+v+2 eq 0 then
        return false, _, _, _;
    end if;

    qtilde := -XF^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*XF
        - (u^2 + u*v + v^2 + u + v + 1);
    f := ((1-u)*XF + 1)*((1-v)*XF + 1)*((u+v+2)*XF + 1)*qtilde;
    yq := u*v*(u+v+1);

    if Degree(f) ne 5 or Discriminant(f) eq 0 or yq eq 0 then
        return false, _, _, _;
    end if;
    if yq^2 ne Evaluate(f, F!-1) then
        return false, _, _, _;
    end if;

    key := Z!u + Characteristic(F)*(Z!v);
    return true, f, yq, key;
end function;

out := Open(output_file, "w");
print "Dumping M_1(8,2,2) halving residue classes";
print "prime_bound", prime_bound, "output_file", output_file;

for p in prime_list do
    F := GF(p);
    PF<XF> := PolynomialRing(F);
    allowed := [];
    nonsingular := 0;
    order8 := 0;

    for u in F do
        for v in F do
            open, f, yq, key := FiniteOpenData(u, v, F, PF);
            if not open then
                continue;
            end if;
            nonsingular +:= 1;
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            DQ := J![XF + F!1, yq];
            if Order(DQ) eq 8 then
                order8 +:= 1;
            end if;
            if IsDivisibleBy2Finite(J, DQ) then
                Append(~allowed, key);
            end if;
        end for;
    end for;

    allowed := Sort(allowed);
    line := IntegerToString(p);
    for key in allowed do
        line cat:= " " cat IntegerToString(key);
    end for;
    Write(out, line cat "\n");
    print "p", p, "open", nonsingular, "Q_order8", order8,
          "allowed_open_halving_residues", #allowed;
end for;

delete out;
quit;
