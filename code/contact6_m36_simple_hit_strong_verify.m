LibraryMode := true;
load "../../verify_simple_torsion_candidate.m";

f66 := 11389248*x^5 - 18252000*x^4 + 42399396*x^3
       - 10288044*x^2 + 29659500*x;
VerifyCandidate("contact6 exact [6,6] seed strong audit", f66, P!0, [[6,6]]);
quit;
