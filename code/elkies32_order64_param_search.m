//////////////////////////////////////////////////////////////////////
//  Targeted search for order 64 on the reconstructed Elkies N=32 family.
//
//  This is not a blind coefficient search.  It searches the genus-0
//  parameter of the reconstructed N=32 component and tests the exact
//  halving equations for D=(0,1)-infinity on residue-sieve survivors.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned B then
    B := 2000;
end if;
if Type(B) eq MonStgElt then
    B := StringToInteger(B);
else
    B := Integers()!B;
end if;
if not assigned output_file then
    output_file := Sprintf("data/elkies32_order64_param_search_B%o.txt", B);
end if;

Q := Rationals();
Zint := Integers();
P2<Z,R,W> := ProjectiveSpace(Q,2);

Fcomp := 3*Z^4*R^4 + 9*Z^4*R^3*W - 16*Z^3*R^4*W
       + 10*Z^4*R^2*W^2 - 56*Z^3*R^3*W^2 + 32*Z^2*R^4*W^2
       + 5*Z^4*R*W^3 - 72*Z^3*R^2*W^3 + 144*Z^2*R^3*W^3
       - 64*Z*R^4*W^3 + Z^4*W^4 - 40*Z^3*R*W^4
       + 208*Z^2*R^2*W^4 - 224*Z*R^3*W^4 + 80*R^4*W^4
       - 8*Z^3*W^5 + 120*Z^2*R*W^5 - 288*Z*R^2*W^5
       + 160*R^3*W^5 + 24*Z^2*W^6 - 160*Z*R*W^6
       + 160*R^2*W^6 - 32*Z*W^7 + 80*R*W^7 + 16*W^8;

C := Curve(P2,Fcomp);
phi := Parametrization(C, C![11/6,1/15,1]);
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

function HasHalfMod(F,z,r,pv)
    if z eq 0 or pv eq 0 then return false; end if;
    a,b,c,pdummy := ABC(z,r);
    f1 := 2*c; f2 := c^2+2*b; f3 := 2*a+2*b*c; f4 := b^2+2*a*c; f5 := a*(2*b-a);
    if f5 eq 0 then return false; end if;
    els := [F!i : i in [0..#F-1]];
    for sgn in [F!-1,F!1] do
        for m in els do
            for n in els do
                al := (f4-m^2)/(2*f5);
                be := (f3-2*m*n-f5*al^2)/(2*f5);
                if f2-n^2-2*m*sgn-2*f5*al*be eq 0 and f1-2*n*sgn-f5*be^2 eq 0 then
                    return true;
                end if;
            end for;
        end for;
    end for;
    return false;
end function;

function HasHalfExact(a,b,c)
    f1 := 2*c; f2 := c^2+2*b; f3 := 2*a+2*b*c; f4 := b^2+2*a*c; f5 := a*(2*b-a);
    if f5 eq 0 then return false, 0, []; end if;
    Rmn<m,n> := PolynomialRing(Q,2);
    for sgn in [-1,1] do
        al := (Rmn!f4-m^2)/(2*f5);
        be := (Rmn!f3-2*m*n-f5*al^2)/(2*f5);
        K2 := Rmn!(f2-n^2-2*m*sgn-2*f5*al*be);
        K1 := Rmn!(f1-2*n*sgn-f5*be^2);
        den := LCM([Denominator(x) : x in Coefficients(K1)] cat [Denominator(x) : x in Coefficients(K2)]);
        I := ideal<Rmn | den*K1, den*K2>;
        if Dimension(I) eq 0 then
            pts := Variety(I);
            if #pts gt 0 then
                return true, sgn, pts;
            end if;
        end if;
    end for;
    return false, 0, [];
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
            den := PDen(zv,rv);
            if den ne 0 then
                pv := 4*PNum(zv,rv)/den;
                if HasHalfMod(Fq,zv,rv,pv) then
                    Append(~allowed, tt);
                end if;
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

out := Open(output_file, "w");
fprintf out, "# Elkies N=32 reconstructed family: search for D divisible by 2\n\n";
fprintf out, "Parameter height bound B = %o\n", B;
fprintf out, "Parametrization coordinate degrees = %o\n", [Degree(e) : e in eqs];
fprintf out, "Sieve primes are applied in the chosen genus-0 parameter.  This is a targeted search, not a proof that no bad-reduction parameter was missed.\n\n";

primes := [11,13,19,23,29,31,37,41,43,47,53,59];
allowed_data := [* *];
fprintf out, "Residue sieve data:\n";
for q in primes do
    A,img := AllowedForPrime(q);
    Append(~allowed_data,<q,A,img>);
    fprintf out, "  p=%o image=%o allowed=%o residues=%o\n", q, img, #A, A;
end for;
fprintf out, "\n";

PX<x> := PolynomialRing(Q);
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
        den := PDen(zv,rv);
        if den eq 0 or zv eq 0 then
            fprintf out, "SURVIVOR t=%o/%o skipped: boundary denominator\n", U,V;
            continue;
        end if;
        pv := 4*PNum(zv,rv)/den;
        if pv eq 0 then
            fprintf out, "SURVIVOR t=%o/%o skipped: p=0\n", U,V;
            continue;
        end if;
        av,bv,cv,pcheck := ABC(zv,rv);
        f := (av*x^3+bv*x^2+cv*x+1)^2 - av^2*x^5*(x+1);
        if Discriminant(f) eq 0 then
            fprintf out, "SURVIVOR t=%o/%o skipped: singular curve\n", U,V;
            continue;
        end if;
        nonsing +:= 1;
        ok,sgn,pts := HasHalfExact(av,bv,cv);
        exact +:= 1;
        fprintf out, "EXACT t=%o/%o result=%o\n", U,V, ok;
        if ok then
            hits +:= 1;
            fprintf out, "  HIT sign=%o points=%o\n", sgn, pts;
            fprintf out, "  z=%o\n  r=%o\n  p=%o\n  a=%o\n  b=%o\n  c=%o\n", zv, rv, pv, av, bv, cv;
        end if;
    end for;
end for;

fprintf out, "\nSUMMARY B=%o total=%o pass_sieve=%o nonsingular=%o exact_checked=%o hits=%o\n", B,total,pass,nonsing,exact,hits;
delete out;
print "Wrote", output_file;
quit;
