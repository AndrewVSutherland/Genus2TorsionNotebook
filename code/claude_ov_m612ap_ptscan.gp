/* claude_ov_m612ap_ptscan.gp   lane 9 ([6,12])   2026-07-25
 * Direct search for rational points of the genus-4 gatekeeper curve
 *    E8 : y^8 + A(x)y^4 + B(x)y^2 + C(x) = 0
 * over every rational x0 = a/b, gcd(a,b)=1, max(|a|,|b|) <= H.
 * claude_ov_m612ap_wfind.m did this to H = 30; here it is pushed to H = 400
 * with a mod-p prefilter (a rational point forces a root of the fibre
 * polynomial mod every good p).  Only x = 0 should survive: the two boundary
 * points b1,b2 of E8 lie over x = 0.
 * NB this GP build does not continue statements across lines in a script
 * file, hence the one-statement-per-line style.
 * Usage: gp -q -f code/claude_ov_m612ap_ptscan.gp
 */
H = 400;
PRS = [7, 11, 13, 17, 23, 29, 31, 37, 41, 43, 47, 53];
Af(X) = 216*X^4 + 72*X^3 - 24*X^2;
Bf(X) = -1296*X^6 - 1728*X^5 - 432*X^4 + 64*X^3;
Cf(X) = -3888*X^8 - 2592*X^7 + 432*X^6 + 288*X^5 - 48*X^4;
fibhit(p, x) = my(a = Mod(Af(x), p), b = Mod(Bf(x), p), c = Mod(Cf(x), p), y2, y4); for (y = 0, p-1, y2 = Mod(y*y, p); y4 = y2*y2; if (y4*y4 + a*y4 + b*y2 + c == 0, return(1))); 0;
rootable = vector(#PRS, i, my(p = PRS[i]); vector(p, k, fibhit(p, k-1)));
for (i = 1, #PRS, print("p = ", PRS[i], " : residues x with a root in the fibre : ", sum(k = 1, PRS[i], rootable[i][k]), "/", PRS[i]));
passes(a, b) = my(p, x); for (i = 1, #PRS, p = PRS[i]; if (b % p == 0, next); x = lift(Mod(a, p) / Mod(b, p)); if (rootable[i][x+1] == 0, return(0))); 1;
cnt = 0; surv = 0; hits = 0;
scan(a, b) = my(x0, f, r); cnt++; if (!passes(a, b), return(0)); surv++; x0 = a/b; f = y^8 + Af(x0)*y^4 + Bf(x0)*y^2 + Cf(x0); r = nfroots(, f); if (#r > 0, hits++; print("HIT x0 = ", x0, "  y = ", r)); 0;
for (b = 1, H, for (a = -H, H, if (gcd(a, b) == 1, scan(a, b))));
print("scanned ", cnt, " coprime pairs with max(|a|,|b|) <= ", H);
print("survived the ", #PRS, "-prime prefilter: ", surv);
print("x0 with a rational point in the fibre: ", hits);
print("PTSCAN_DONE");
quit;
