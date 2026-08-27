// Lane 8 (overnight 2026-07-25), part B: recover explicit genus-2 models for the
// remaining HLP split targets that the repository has never written down --
//   Z/5 x Z/10  (HLP sec 3.5, N=N'=10:  Delta_10(t) y^2 = Delta_10(u))
//   Z/2 x Z/24  (HLP sec 3.4, E_{2,6}^t glued to E_{2,8}^u; 2-torsion trivially matched)
//   Z/3 x Z/12  (HLP sec 3.5, N=12, N'=6:  Delta_12(t) y^2 = (u+1)(9u+1), a conic)
// and to build FAMILIES (many members, not one seed) for
//   Z/35        (E_7^{-1} glued to E_5^u,  u from the rank-2 curve y^2 = -26 u(u^2-11u-1))
//   Z/45        (E_9^{-5} glued to E_5^u,  u from the rank-2 curve y^2 = -930 u(u^2-11u-1))
//
// Uses Fork() to fan the per-fibre Prop-4 + TorsionSubgroup work over NCH children.
SetColumns(0);
if not assigned NCH then NCH := 16; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

QQ := Rationals();
P<x> := PolynomialRing(QQ);

function BC(N, t)
    case N:
        when 5:  return t, t;
        when 6:  return t^2+t, t;
        when 7:  return t^3-t^2, t^2-t;
        when 9:  return t^5-2*t^4+2*t^3-t^2, t^3-t^2;
        when 10: return (2*t^5-3*t^4+t^3)/(t^2-3*t+1)^2, (-2*t^3+3*t^2-t)/(t^2-3*t+1);
        when 12: return (12*t^6-30*t^5+34*t^4-21*t^3+7*t^2-t)/(t-1)^4,
                        (-6*t^4+9*t^3-5*t^2+t)/(t-1)^3;
        when 26: return (-2*t^3+14*t^2-22*t+10)/((t+3)^2*(t-3)^2), (-2*t+10)/((t+3)*(t-3));
        when 28: return (16*t^3+16*t^2+6*t+1)/(8*t^2-1)^2,
                        (16*t^3+16*t^2+6*t+1)/(2*t*(4*t+1)*(8*t^2-1));
    end case;
end function;

function UnivCubic(N, t)
    b, c := BC(N, t);
    return x^3 - b*x^2 + ((1-c)*x - b)^2/4;
end function;

function Prop4(f, g)
    L := SplittingField(f*g);
    PL<X> := PolynomialRing(L);
    ra := [r[1] : r in Roots(PL!f)];
    rb := [r[1] : r in Roots(PL!g)];
    if #ra ne 3 or #rb ne 3 then return []; end if;
    Df := QQ!Discriminant(f);  Dg := QQ!Discriminant(g);
    out := [];
    for perm in Permutations({1,2,3}) do
        al := ra;
        be := [rb[perm[1]], rb[perm[2]], rb[perm[3]]];
        a1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1])
            + (al[1]-al[3])^2/(be[1]-be[3]);
        b1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1])
            + (be[1]-be[3])^2/(al[1]-al[3]);
        a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
        b2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
        if a1 eq 0 or a2 eq 0 or b1 eq 0 or b2 eq 0 then continue; end if;
        A := Dg*a1/a2;  B := Df*b1/b2;
        hh := -( A*(al[2]-al[1])*(al[1]-al[3])*X^2 + B*(be[2]-be[1])*(be[1]-be[3]))
            * ( A*(al[3]-al[2])*(al[2]-al[1])*X^2 + B*(be[3]-be[2])*(be[2]-be[1]))
            * ( A*(al[1]-al[3])*(al[3]-al[2])*X^2 + B*(be[1]-be[3])*(be[3]-be[2]));
        ok := true; csQ := [];
        for cc in Coefficients(hh) do
            fl, q := IsCoercible(QQ, cc);
            if not fl then ok := false; break; end if;
            Append(~csQ, q);
        end for;
        if not ok then continue; end if;
        hQ := P!csQ;
        if Degree(hQ) ne 6 or not IsSquarefree(hQ) then continue; end if;
        Append(~out, hQ);
    end for;
    return out;
end function;

function CleanSextic(h)
    d := LCM([Denominator(c) : c in Coefficients(h)]);
    h2 := P!(d^2*h);
    cont := GCD([Integers()!c : c in Coefficients(h2)]);
    sq := 1;
    for pr in Factorization(cont) do sq *:= pr[1]^(2*(pr[2] div 2)); end for;
    return P!(h2/sq);
end function;

// -------------------------------------------------------------------------
// build the job list  <tag, f, g>  in the PARENT (cheap; children inherit it)
// -------------------------------------------------------------------------
jobs := [];

// ---- Z/3 x Z/12 : Delta_12(t) y^2 = (u+1)(9u+1).  Conic through (u,y)=(-1,0):
//      u = -1 + s,  s = 8/(9 - Delta_12(t) m^2),  y = m s.
D12 := func<t | (2*t^2-2*t+1)*(6*t^2-6*t+1)>;
for t in [QQ!2, QQ!3, QQ!(1/2), QQ!(-1), QQ!(1/3), QQ!(3/2), QQ!4, QQ!(-2)] do
    for m in [QQ!1, QQ!2, QQ!(1/2), QQ!3, QQ!(1/3), QQ!(-1), QQ!(-2), QQ!(2/3)] do
        den := 9 - D12(t)*m^2;
        if den eq 0 then continue; end if;
        u := -1 + 8/den;
        if u in {QQ!0, QQ!(-1), QQ!(-1/9)} then continue; end if;
        f := UnivCubic(12, t); g := UnivCubic(6, u);
        if Discriminant(f) eq 0 or Discriminant(g) eq 0 then continue; end if;
        Append(~jobs, <Sprintf("Z3xZ12 t=%o m=%o u=%o", t, m, u), f, g>);
    end for;
end for;

// ---- Z/2 x Z/24 : E_{2,6}^t glued to E_{2,8}^u -- both have full rational 2-torsion,
//      so EVERY pair works (HLP sec 3.4); the group depends on the labelling.
for t in [QQ!1, QQ!2, QQ!(1/2), QQ!5, QQ!(-1), QQ!(1/3), QQ!7, QQ!(2/3)] do
    for u in [QQ!1, QQ!2, QQ!(1/2), QQ!(1/3), QQ!3, QQ!(-1), QQ!(1/4), QQ!(3/2)] do
        f := UnivCubic(26, t); g := UnivCubic(28, u);
        if Discriminant(f) eq 0 or Discriminant(g) eq 0 then continue; end if;
        if jInvariant(EllipticCurve(f)) eq jInvariant(EllipticCurve(g)) then continue; end if;
        Append(~jobs, <Sprintf("Z2xZ24 t=%o u=%o", t, u), f, g>);
    end for;
end for;

// ---- Z/5 x Z/10 : Delta_10(t) y^2 = Delta_10(u), Delta_10(v) = 8v^3-8v^2+1.
//      For each t the (u,y)-curve is genus 1 with the rational point (t,1); take
//      multiples of that point on the elliptic curve to obtain u != t.
D10 := func<v | 8*v^3 - 8*v^2 + 1>;
for t in [QQ!2, QQ!3, QQ!(-1), QQ!(1/3), QQ!4, QQ!(3/2), QQ!5, QQ!(-2), QQ!(1/2), QQ!(-1/2)] do
    Dt := D10(t);
    if Dt eq 0 then continue; end if;
    // Y^2 = Dt*(8u^3-8u^2+1) with Y = Dt*y ; rational point (u,Y) = (t, Dt).
    // Direct point search (NEVER the birational-map preimage: Magma's @@ on a
    // CurveQuotient/EllipticCurve map picks up patch-indeterminacy junk).
    dsc := Numerator(Dt)*Denominator(Dt);      // same square class as Dt, integral
    Cg1 := HyperellipticCurve(dsc*(8*x^3-8*x^2+1));
    pts := Points(Cg1 : Bound := 400);
    seen := {t};
    for pp in pts do
        if pp[3] eq 0 then continue; end if;
        uu := pp[1]/pp[3];
        if uu in seen then continue; end if;
        Include(~seen, uu);
        if D10(uu) eq 0 then continue; end if;
        f := UnivCubic(10, t); g := UnivCubic(10, uu);
        if Discriminant(f) eq 0 or Discriminant(g) eq 0 then continue; end if;
        if jInvariant(EllipticCurve(f)) eq jInvariant(EllipticCurve(g)) then continue; end if;
        Append(~jobs, <Sprintf("Z5xZ10 t=%o u=%o", t, uu), f, g>);
    end for;
end for;

// ---- Z/35 family : E_7^{-1} glued to E_5^u for u on  y^2 = -26 u(u^2-11u-1)
Cq := HyperellipticCurve(-26*(x^3 - 11*x^2 - x));
us35 := {QQ!(1/26)};
for pp in Points(Cq : Bound := 2000) do
    if pp[3] eq 0 then continue; end if;
    uu := pp[1]/pp[3];
    if uu eq 0 then continue; end if;
    Include(~us35, uu);
end for;
printf "Z35 candidate u values: %o\n", #us35;
f7m1 := UnivCubic(7, QQ!(-1));
for uu in us35 do
    g := UnivCubic(5, uu);
    if Discriminant(g) eq 0 then continue; end if;
    Append(~jobs, <Sprintf("Z35 u=%o", uu), f7m1, g>);
end for;

// ---- Z/45 family : E_9^{-5} glued to E_5^u for u on  y^2 = -930 u(u^2-11u-1)
Cq45 := HyperellipticCurve(-930*(x^3 - 11*x^2 - x));
us45 := {QQ!(93/10)};
for pp in Points(Cq45 : Bound := 2000) do
    if pp[3] eq 0 then continue; end if;
    uu := pp[1]/pp[3];
    if uu eq 0 then continue; end if;
    Include(~us45, uu);
end for;
printf "Z45 candidate u values: %o\n", #us45;
f9m5 := UnivCubic(9, QQ!(-5));
for uu in us45 do
    g := UnivCubic(5, uu);
    if Discriminant(g) eq 0 then continue; end if;
    Append(~jobs, <Sprintf("Z45 u=%o", uu), f9m5, g>);
end for;

printf "JOBS %o\n", #jobs;

dir := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/lane8b";
System("rm -rf " cat dir cat "; mkdir -p " cat dir);

for c in [0..NCH-1] do
    pid := Fork();
    if pid eq 0 then
        Fh := Open(dir cat "/part" cat IntegerToString(c) cat ".txt", "w");
        for i in [c+1..#jobs by NCH] do
            tag := jobs[i][1]; f := jobs[i][2]; g := jobs[i][3];
            ok := true; hs := [];
            try hs := Prop4(f, g); catch e ok := false; end try;
            if not ok then Puts(Fh, tag cat " | ERROR"); continue; end if;
            if #hs eq 0 then Puts(Fh, tag cat " | NO-RATIONAL-GLUING"); continue; end if;
            for h in hs do
                hc := CleanSextic(h);
                try
                    C := HyperellipticCurve(hc);
                    J := Jacobian(C);
                    T := TorsionSubgroup(J);
                    inv := Invariants(T);
                    E := EllipticCurve(f); F := EllipticCurve(g);
                    Puts(Fh, tag cat " | TOR " cat Sprint(inv) cat " | ORD " cat IntegerToString(&*[Integers()|1] * (#T))
                             cat " | 2RK " cat Sprint(Invariants(TwoTorsionSubgroup(J)))
                             cat " | E " cat IntegerToString(Conductor(E)) cat " " cat Sprint(Invariants(TorsionSubgroup(E)))
                             cat " | F " cat IntegerToString(Conductor(F)) cat " " cat Sprint(Invariants(TorsionSubgroup(F)))
                             cat " | MODEL " cat Sprint(hc));
                catch e
                    Puts(Fh, tag cat " | TORSION-ERROR | MODEL " cat Sprint(hc));
                end try;
            end for;
        end for;
        Flush(Fh); delete Fh;
        quit;
    end if;
end for;
WaitForAllChildren();

for c in [0..NCH-1] do
    fn := dir cat "/part" cat IntegerToString(c) cat ".txt";
    for l in Split(Read(fn), "\n") do
        if #l gt 0 then print l; end if;
    end for;
end for;
print "SEARCH_DONE prop4b";
quit;
