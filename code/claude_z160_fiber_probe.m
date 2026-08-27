// claude_z160_fiber_probe.m — understand the UNEXPECTED dim-1 component of the
// second-5-contact system on the Elkies [32] family, by finite-field fibers:
// for u0 in F_p, solve the 0-dim system {ftil = h5^2 - c5 (x-rho)^5 coefficients}
// in (h2,h1,h0,rho) over F_p and inspect the solutions' shape.  Also gate-count
// 5 | #J(F_p) as the necessary condition for a rational 5-torsion class.
SetColumns(0);
SetMemoryLimit(4*10^9);
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
Q := Rationals();
cfs := [ Numerator(Coefficient(ftil, i)) : i in [0..5] ];

for p in [101, 211] do
  Fp := GF(p);
  Pp<t> := PolynomialRing(Fp);
  cfp := [ Pp!c : c in cfs ];
  R4<h2, h1, h0, rho> := PolynomialRing(Fp, 4);
  Pxx<X> := PolynomialRing(R4);
  nsol := 0; nfib := 0; ngate := 0;
  shapes := [];
  for u0 in Fp do
    fv := [ Evaluate(c, u0) : c in cfp ];
    if fv[6] eq 0 then continue; end if;     // degenerate member
    fpol := Pp![ fv[i] : i in [1..6] ];
    if Discriminant(fpol) eq 0 then continue; end if;
    nfib +:= 1;
    // gate: 5 | #J(F_p)
    np := #Jacobian(HyperellipticCurve(fpol));
    if np mod 5 eq 0 then ngate +:= 1; end if;
    // contact fiber
    Fx := &+[ Pxx | fv[i+1]*X^i : i in [0..5] ];
    H5 := h2*X^2 + h1*X + h0;
    c5 := -fv[6];
    Res := Fx - H5^2 + c5*(X - rho)^5;
    I := ideal< R4 | [Coefficient(Res, i) : i in [0..4]] >;
    if Dimension(I) ge 0 then
      V := Variety(I);
      if #V gt 0 then
        nsol +:= 1;
        if #shapes lt 6 then
          for v in V do
            rr := v[4];
            info := Sprintf("u0=%o sol (h2,h1,h0,rho)=%o f(rho)=%o h5(rho)^2=%o",
              u0, v, Evaluate(fpol, rr), (v[1]*rr^2+v[2]*rr+v[3])^2);
            Append(~shapes, info);
          end for;
        end if;
      end if;
    end if;
  end for;
  printf "p=%o: fibers=%o, contact-solvable=%o, gate 5|#J=%o\n", p, nfib, nsol, ngate;
  for s in shapes do printf "  %o\n", s; end for;
end for;
print "ALL_DONE";
quit;
