//////////////////////////////////////////////////////////////////////
//  Search for [2,32] on the reconstructed Elkies N=32 family.
//
//  Algebraic condition: the reconstructed family already has one rational
//  finite Weierstrass root r with (r,0)-oo = 16*((0,1)-oo).  To get an
//  independent rational 2-torsion class, force a second rational root s != r
//  of f(x).  This script searches the genus-0 parameter of the reconstructed
//  N=32 component, using finite-field extra-root residues before exact
//  factorization.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned B then
    B := 5000;
end if;
if Type(B) eq MonStgElt then
    B := StringToInteger(B);
else
    B := Integers()!B;
end if;
if not assigned max_hits then
    max_hits := 20;
end if;
if Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned output_file then
    output_file := Sprintf("data/elkies32_extra2_param_search_B%o.txt", B);
end if;

Q := Rationals();
Zint := Integers();
P<x> := PolynomialRing(Q);
P2<Z,R,W> := ProjectiveSpace(Q,2);

Fcomp := 3*Z^4*R^4 + 9*Z^4*R^3*W - 16*Z^3*R^4*W
       + 10*Z^4*R^2*W^2 - 56*Z^3*R^3*W^2 + 32*Z^2*R^4*W^2
       + 5*Z^4*R*W^3 - 72*Z^3*R^2*W^3 + 144*Z^2*R^3*W^3
       - 64*Z*R^4*W^3 + Z^4*W^4 - 40*Z^3*R*W^4
       + 208*Z^2*R^2*W^4 - 224*Z*R^3*W^4 + 80*R^4*W^4
       - 8*Z^3*W^5 + 120*Z^2*R*W^5 - 288*Z*R^2*W^5
       + 160*R^3*W^5 + 24*Z^2*W^6 - 160*Z*R*W^6
       + 160*R^2*W^6 - 32*Z*W^7 + 80*R*W^7 + 16*W^8;

Cbase := Curve(P2,Fcomp);
phi := Parametrization(Cbase, Cbase![11/6,1/15,1]);
eqs := DefiningEquations(phi);

function PNum(z,r)
    return (z - 2)*(r + 1)^2*(z^2*r + z^2 - 4*r)
         * (z^4*r^2 + 2*z^4*r - 10*z^3*r^2 + z^4 - 16*z^3*r
            + 28*z^2*r^2 - 6*z^3 + 44*z^2*r - 16*z*r^2
            + 16*z^2 - 56*z*r + 16*r^2 - 24*z + 32*r + 16);
end function;

function PDen(z,r)
    return r*z*(z^8*r^4 + 4*z^8*r^3 - 24*z^7*r^4 + 48*z^6*r^5
        + 6*z^8*r^2 - 80*z^7*r^3 + 284*z^6*r^4 - 256*z^5*r^5
        + 4*z^8*r - 100*z^7*r^2 + 632*z^6*r^3 - 1200*z^5*r^4
        + 320*z^4*r^5 + z^8 - 56*z^7*r + 672*z^6*r^2
        - 2432*z^5*r^3 + 2224*z^4*r^4 - 12*z^7 + 344*z^6*r
        - 2496*z^5*r^2 + 5312*z^4*r^3 - 2304*z^3*r^4 + 68*z^6
        - 1248*z^5*r + 5808*z^4*r^2 - 7232*z^3*r^3 + 1664*z^2*r^4
        - 240*z^5 + 2976*z^4*r - 8832*z^3*r^2 + 6272*z^2*r^3
        - 768*z*r^4 + 576*z^4 - 4800*z^3*r + 8640*z^2*r^2
        - 3328*z*r^3 + 256*r^4 - 960*z^3 + 5120*z^2*r
        - 5120*z*r^2 + 1024*r^3 + 1088*z^2 - 3328*z*r
        + 1536*r^2 - 768*z + 1024*r + 256);
end function;

function ABC(zv, rv)
    pv := 4*PNum(zv,rv)/PDen(zv,rv);
    av := zv*pv;
    cv := (-pv*zv^4 + 4*pv*zv^3 - 8*pv*zv^2 + 8*zv^2 + 8*pv*zv - 16)/(4*zv^2);
    bv := pv*(zv-1) + cv - 1;
    return av,bv,cv,pv;
end function;

function FPolynomial(a,b,c, Rpoly)
    X := Rpoly.1;
    return (a*X^3 + b*X^2 + c*X + 1)^2 - a^2*X^5*(X+1);
end function;

function ExtraRootMod(Fq, z, r)
    if z eq 0 then return false; end if;
    den := PDen(z,r);
    if den eq 0 then return false; end if;
    pval := 4*PNum(z,r)/den;
    if pval eq 0 then return false; end if;
    a,b,c,pdummy := ABC(z,r);
    Rpoly<X> := PolynomialRing(Fq);
    f := FPolynomial(a,b,c,Rpoly);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then return false; end if;
    for s in Fq do
        if s ne r and Evaluate(f,s) eq 0 then
            return true;
        end if;
    end for;
    return false;
end function;

function AllowedForPrime(q)
    Fq := GF(q);
    PR := PolynomialRing(Fq,2);
    eqsf := [PR!e : e in eqs];
    allowed := [];
    img := {};
    for tt in [0..q-1] cat [-1] do
        uu := tt eq -1 select Fq!1 else Fq!tt;
        vv := tt eq -1 select Fq!0 else Fq!1;
        Zv := Evaluate(eqsf[1],[uu,vv]);
        Rv := Evaluate(eqsf[2],[uu,vv]);
        Wv := Evaluate(eqsf[3],[uu,vv]);
        if Wv ne 0 then
            zv := Zv/Wv;
            rv := Rv/Wv;
            Include(~img,<zv,rv>);
            if ExtraRootMod(Fq,zv,rv) then
                Append(~allowed, tt);
            end if;
        end if;
    end for;
    return allowed, #img;
end function;

function PassSieve(U,V, allowed_data)
    for rec in allowed_data do
        q := rec[1];
        A := rec[2];
        if (V mod q) eq 0 then
            tt := -1;
        else
            tt := Zint!((GF(q)!U)/(GF(q)!V));
        end if;
        if not tt in A then
            return false;
        end if;
    end for;
    return true;
end function;

function ExtraRootsExact(f, r)
    roots := [];
    for pair in Factorization(f) do
        g := pair[1];
        if Degree(g) eq 1 then
            s := -Coefficient(g,0)/Coefficient(g,1);
            if s ne r then
                Append(~roots, s);
            end if;
        end if;
    end for;
    return roots;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, P!0;
end function;

function VerifyTorsion(f, r, extra_roots)
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D := J![x, 1];
    ordD := Order(D);
    Tbuilt := J![x-r, 0];
    ok32 := ordD eq 32 and 16*D eq Tbuilt;
    independent := false;
    good_s := [];
    zero := J!0;
    for s in extra_roots do
        T2 := J![x-s, 0];
        if T2 ne zero and T2 ne Tbuilt then
            independent := true;
            Append(~good_s, s);
        end if;
    end for;
    return ok32, independent, ordD, good_s;
end function;

out := Open(output_file, "w");
fprintf out, "# Elkies N=32 reconstructed family: search for extra independent 2-torsion ([2,32])\n\n";
fprintf out, "Parameter height bound B = %o\n", B;
fprintf out, "Parametrization coordinate degrees = %o\n", [Degree(e) : e in eqs];
fprintf out, "Algebraic condition: f(s)=0 for a rational s distinct from the built-in root r.\n\n";

diagnostic_primes := [11,13];
primes := [19,23,29,31,37,41,43,47,53,59,61,67,71,73];
allowed_data := [* *];
fprintf out, "Diagnostic primes not used for sieve:\n";
for q in diagnostic_primes do
    A,img := AllowedForPrime(q);
    fprintf out, "  p=%o image=%o allowed=%o residues=%o\n", q, img, #A, A;
end for;
fprintf out, "Residue sieve data:\n";
for q in primes do
    A,img := AllowedForPrime(q);
    Append(~allowed_data,<q,A,img>);
    fprintf out, "  p=%o image=%o allowed=%o residues=%o\n", q, img, #A, A;
end for;
fprintf out, "\n";

total := 0; pass := 0; nonsing := 0; exact := 0; hits := 0;

for V in [1..B] do
    for U in [-B..B] do
        if GCD(U,V) ne 1 then continue; end if;
        total +:= 1;
        if not PassSieve(U,V,allowed_data) then continue; end if;
        pass +:= 1;
        Zv := Evaluate(eqs[1],[Q!U,Q!V]);
        Rv := Evaluate(eqs[2],[Q!U,Q!V]);
        Wv := Evaluate(eqs[3],[Q!U,Q!V]);
        if Wv eq 0 then
            fprintf out, "SURVIVOR t=%o/%o skipped: W=0\n", U,V;
            continue;
        end if;
        zv := Zv/Wv;
        rv := Rv/Wv;
        if zv eq 0 or PDen(zv,rv) eq 0 then
            fprintf out, "SURVIVOR t=%o/%o skipped: boundary denominator\n", U,V;
            continue;
        end if;
        av,bv,cv,pv := ABC(zv,rv);
        if av eq 0 then
            fprintf out, "SURVIVOR t=%o/%o skipped: a=0\n", U,V;
            continue;
        end if;
        f := FPolynomial(av,bv,cv,P);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            fprintf out, "SURVIVOR t=%o/%o skipped: singular/bad degree\n", U,V;
            continue;
        end if;
        nonsing +:= 1;
        roots := ExtraRootsExact(f,rv);
        exact +:= 1;
        fprintf out, "EXACT t=%o/%o extra_roots=%o\n", U,V, roots;
        if #roots gt 0 then
            ok32, independent, ordD, good_s := VerifyTorsion(f,rv,roots);
            simple, pcert, Lp := SimpleCertificate(f);
            fprintf out, "  VERIFY ok32=%o independent=%o ordD=%o good_extra_roots=%o simple_cert=%o pcert=%o Lp=%o\n", ok32, independent, ordD, good_s, simple, pcert, Lp;
            fprintf out, "  z=%o\n  r=%o\n  p=%o\n  a=%o\n  b=%o\n  c=%o\n  f=%o\n", zv,rv,pv,av,bv,cv,f;
            if ok32 and independent then
                hits +:= 1;
            end if;
            if hits ge max_hits then
                fprintf out, "Reached max_hits=%o\n", max_hits;
                fprintf out, "\nSUMMARY B=%o total=%o pass_sieve=%o nonsingular=%o exact_checked=%o hits=%o\n", B,total,pass,nonsing,exact,hits;
                delete out;
                print "Wrote", output_file;
                quit;
            end if;
        end if;
    end for;
end for;

fprintf out, "\nSUMMARY B=%o total=%o pass_sieve=%o nonsingular=%o exact_checked=%o hits=%o\n", B,total,pass,nonsing,exact,hits;
delete out;
print "Wrote", output_file;
quit;
