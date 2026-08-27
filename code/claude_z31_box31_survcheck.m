// Exact Magma stage for claude_z31_box31 survivors (Task B5 validation).
// Usage: magma -b C6:=-38 C5:=31 C4:=-6 C3:=-29 C2:=-40 C1:=9 C0:=2 claude_z31_box31_survcheck.m
// Defaults are the first survivor found by the spot-11 H=40 bench shard.
SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned C6 then C6 := -38; elif Type(C6) eq MonStgElt then C6 := StringToInteger(C6); end if;
if not assigned C5 then C5 := 31;  elif Type(C5) eq MonStgElt then C5 := StringToInteger(C5); end if;
if not assigned C4 then C4 := -6;  elif Type(C4) eq MonStgElt then C4 := StringToInteger(C4); end if;
if not assigned C3 then C3 := -29; elif Type(C3) eq MonStgElt then C3 := StringToInteger(C3); end if;
if not assigned C2 then C2 := -40; elif Type(C2) eq MonStgElt then C2 := StringToInteger(C2); end if;
if not assigned C1 then C1 := 9;   elif Type(C1) eq MonStgElt then C1 := StringToInteger(C1); end if;
if not assigned C0 then C0 := 2;   elif Type(C0) eq MonStgElt then C0 := StringToInteger(C0); end if;
P<x> := PolynomialRing(Rationals());
f := C6*x^6 + C5*x^5 + C4*x^4 + C3*x^3 + C2*x^2 + C1*x + C0;
printf "SURVCHECK f = %o\n", f;
C := HyperellipticCurve(f);
J := Jacobian(C);
t0 := Cputime();
T := TorsionSubgroup(J);
inv := Invariants(T);
printf "SURVCHECK torsion invariants = %o  (%o s)\n", inv, Cputime(t0);
printf "SURVCHECK has31 = %o\n", exists{d : d in inv | d mod 31 eq 0};
quit;
