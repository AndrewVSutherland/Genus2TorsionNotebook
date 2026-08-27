// Split-census leg: systematic bielliptic genus-2 curves
//   C: y^2 = a*x^6 + b*x^4 + c*x^2 + d   (double cover of two elliptic curves)
// Jacobian splits by construction (up to isogeny, over Q).  Funnel:
// gcd(#J(F_p)) prefilter at p in {5,7,11} (threshold), then exact torsion.
// Params: H (box), Sh (shard 0..NSh-1 on coefficient a), NSh.
SetColumns(0);
if not assigned H then H := 12; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned Sh then Sh := 0; elif Type(Sh) eq MonStgElt then Sh := StringToInteger(Sh); end if;
if not assigned NSh then NSh := 3; elif Type(NSh) eq MonStgElt then NSh := StringToInteger(NSh); end if;
if not assigned MinOrd then MinOrd := 12; elif Type(MinOrd) eq MonStgElt then MinOrd := StringToInteger(MinOrd); end if;
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
P<x> := PolynomialRing(Rationals());
seen := {};
nproc := 0; nex := 0;
for a in [-H..H] do
    if a eq 0 then continue; end if;
    if ((a + H) mod NSh) ne Sh then continue; end if;
    for b in [-H..H], c in [-H..H], d in [1..H] do
        // d>0 WLOG up to twist y->iy? no; keep d in 1..H plus separately d in -H..-1 
        // handled via (a,b,c,d) -> (-a,-b,-c,-d) twist NOT torsion-equal; keep both signs:
        for dd in [d, -d] do
            f := a*x^6 + b*x^4 + c*x^2 + dd;
            if not IsSquarefree(f) then continue; end if;
            nproc +:= 1;
            // prefilter
            g := 0;
            ok := true;
            for p in [5,7,11,13] do
                if Integers()!Discriminant(f) mod p eq 0 then continue; end if;
                Cp := HyperellipticCurve(PolynomialRing(GF(p))!f);
                np := #Jacobian(Cp);
                g := GCD(g, np);
                if g lt MinOrd then ok := false; break; end if;
            end for;
            if not ok or g lt MinOrd then continue; end if;
            C := HyperellipticCurve(f);
            ih := Hash(G2Invariants(C));
            if ih in seen then continue; end if;
            Include(~seen, ih);
            nex +:= 1;
            T := TorsionSubgroup(Jacobian(C));
            if #T ge MinOrd then
                printf "HIT [%o,%o,%o,%o] TORSION %o order %o\n", a,b,c,dd, Invariants(T), #T;
            end if;
        end for;
    end for;
    printf "PROGRESS a=%o processed %o exact %o\n", a, nproc, nex;
end for;
printf "SEARCH_DONE shard %o/%o H=%o processed %o exact %o\n", Sh, NSh, H, nproc, nex;
quit;
