// Lane 2: scan the fibration of the contact-7 THREE-ROOT surface S by z_1 = c.
// Generic fibre is a (3,3) curve of genus 4; c=2 degenerates to genus 1 (rank-1 elliptic).
// Find every other c of small height whose fibre has genus <= 1  (=> potential infinite family).
SetColumns(0);
Q := Rationals();
R<x,y> := PolynomialRing(Q,2);

// Degenerate z-values (never a genuine rational Weierstrass point of the chart quintic):
//   z=1  (v=0,  r=1 is a double root),  z=1/2 (v=1, r=0, divided out by the x^2),
//   z=0/oo (v=oo/-1).  c=1/2 also makes T=1, which degenerates the pencil.
function FibreCurve(cv)
  if cv eq 1 or cv eq 1/2 or cv eq 0 then return 0, 0, 0; end if;
  T := 3 - 1/(1-cv);
  if T eq 1 then return 0, 0, 0; end if;
  e1 := 6 - cv; e4 := 2/cv;
  A := (T-2)/(T-1);
  B := (T*(cv-5+2/cv) - 3*cv + 14)/(T-1);
  Nm := func<t| -(t^4 - e1*t^3 - B*t + e4)>;
  Dn := func<t| t^2 - A*t>;
  N := Nm(x)*Dn(y) - Nm(y)*Dn(x);
  if N eq 0 then return 0, A, B; end if;
  F := ExactQuotient(N, x-y);
  return F, A, B;
end function;

cands := [];
for d in [1..8] do
  for n in [-40..40] do
    if n eq 0 then continue; end if;
    if GCD(Abs(n),d) ne 1 then continue; end if;
    cv := Q!n/d;
    if cv eq 1 or cv eq 1/2 then continue; end if;   // degenerate z-values
    Append(~cands, cv);
  end for;
end for;
// plus every z-value occurring in the eleven recorded three-root hits
extra := [ Q| -1/2,4,5/2,-1/9,2,-7/3,-1/4,-8/7,22/7,19/4,9/5,17/21,25/21,18/13,49/39,
            -61/450,625/114,-41/124,-16/17,289/124,297/133,361/525,50/33,189/155,121/155,49/36,50/63 ];
for e in extra do if not (e in cands) then Append(~cands,e); end if; end for;
printf "scanning %o values of c\n", #cands; 

low := [];
for cv in cands do
  F, A, B := FibreCurve(cv);
  if F cmpeq 0 then printf "c=%o : degenerate (F=0)\n", cv; continue; end if;
  fa := Factorisation(F);
  tot := 0; gl := [];
  for g in fa do
    C := Curve(AffineSpace(R), g[1]);
    gg := -1;
    try gg := Genus(ProjectiveClosure(C)); catch e gg := -1; end try;
    Append(~gl, <Degree(g[1],x), Degree(g[1],y), g[2], gg>);
  end for;
  mg := Min([ t[4] : t in gl ]);
  if #fa gt 1 or mg le 1 then
    printf "c=%o  A=%o B=%o  factors=%o  LOWGENUS\n", cv, A, B, gl;
    Append(~low, <cv, gl>);
  end if;
end for;
printf "--- %o low-genus / reducible fibres found ---\n", #low;
for t in low do printf "  c=%o  %o\n", t[1], t[2]; end for;
printf "CSCAN_DONE\n";
quit;
