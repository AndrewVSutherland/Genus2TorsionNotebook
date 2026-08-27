//////////////////////////////////////////////////////////////////////
// E8 (m612 relative-3 gatekeeper) in Tuitman's Coleman machinery.
// Stage 1: coleman_data at a good prime; rational point search;
// parity split of the basis under the Prym involution y -> -y.
//////////////////////////////////////////////////////////////////////
load "coleman.m";
Q := y^8 + (216*x^4+72*x^3-24*x^2)*y^4 + (-1296*x^6-1728*x^5-432*x^4+64*x^3)*y^2
     + (-3888*x^8-2592*x^7+432*x^6+288*x^5-48*x^4);
print "model set: E8, monic degree 8 in y, even in y";
for p in [7,13,31,37] do
    printf "=== trying p = %o ===\n", p;
    try
        tt := Cputime();
        data := coleman_data(Q, p, 15);
        printf "coleman_data OK at p=%o (%.1os): genus g = %o\n", p, Cputime(tt), data`g;
        // rational points
        qpts := Q_points(data, 10000);
        printf "Q_points (bound 1e4): %o found\n", #qpts;
        for P in qpts do
            printf "  x=%o inf=%o b1=%o\n", P`x, P`inf, P`b[1];
        end for;
        // basis parity under y -> -y : examine data`basis entries
        print "basis vectors (coeffs w.r.t. integral basis):";
        for i in [1..#data`basis] do
            printf "  basis[%o] = %o\n", i, data`basis[i];
        end for;
        print "W0 ="; print data`W0;
        break;
    catch e
        printf "p=%o failed: %o\n", p, e`Object;
    end try;
end for;
print "E8_COLEMAN1_DONE";
quit;
