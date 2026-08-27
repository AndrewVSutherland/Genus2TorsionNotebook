/* Lane 8 (2026-07-25, resumed session): CORRECTED + EXTENDED landscape census.
 *
 * Fixes two defects of code/claude_ov_lane8_landscape.gp:
 *   (1) it tested only  N | #J(F_p).  For the targets whose group is NOT cyclic
 *       ([7,7], [5,10], [2,24], [3,12]) that is the wrong condition:
 *       (Z/n)^2 <= J(F_p) needs n-RANK >= 2, i.e. (T-1)^2 | chi(T) mod n
 *       together with n^2 | chi(1).
 *   (2) it only ran at p = 3,5,7,11,13 (full curve enumeration, O(p^3) work).
 *       Here we enumerate ISOGENY CLASSES -- all Weil-admissible (s1,s2) --
 *       which is O(p) classes per prime, so p can reach several hundred.
 *
 * Weil-admissibility is tested exactly, in integers:
 *   chi = (T^2-aT+p)(T^2-bT+p) with a+b=s1, ab=s2-2p, a,b real in [-2sqrt p,2sqrt p]
 *   <=>  s1^2-4*s2+8*p >= 0,  s2+2*p >= 0,  (s2+2*p)^2 >= 4*p*s1^2,  s1^2 <= 16*p.
 * This is a SUPERSET of the isogeny classes of abelian surfaces (a fortiori of
 * Jacobians), so "no compatible class at p" is SAFE (forced bad reduction),
 * while "a compatible abs-simple class exists" is an upper bound.
 */

abssimple(chi) = {my(v); if(!polisirreducible(chi), return(0)); for(n=2,12, v = polresultant(subst(chi,x,y), x-y^n, y); v = v/polcoeff(v,poldegree(v)); if(poldegree(v)!=4 || !polisirreducible(v), return(0))); 1;}

targets = [["Z/35    ",35,[]], ["Z/45    ",45,[]], ["Z/63    ",63,[]], ["Z/70    ",70,[]], ["Z/7xZ/7 ",49,[7]], ["Z/5xZ/10",50,[5]], ["Z/2xZ/24",48,[2]], ["Z/3xZ/12",36,[3]], ["Z/40 CTL",40,[]], ["Z/30 CTL",30,[]], ["Z/48 CTL",48,[]]];

rank2ok(chi,n,nJ) = {if(nJ % (n^2) != 0, return(0)); ((Mod(1,n)*chi) % (Mod(1,n)*(x-1)^2)) == 0;}

weilok(p,s1,s2) = {(s1^2 <= 16*p) && (s1^2 - 4*s2 + 8*p >= 0) && (s2 + 2*p >= 0) && ((s2+2*p)^2 >= 4*p*s1^2);}

census(plist) = {my(certprimes = vector(#targets,i,[]), badprimes = vector(#targets,i,[])); for(pi=1,#plist, my(p=plist[pi], s1max=sqrtint(16*p), ncls=0, comp=vector(#targets), simp=vector(#targets), ex=vector(#targets,i,[])); for(s1=-s1max,s1max, for(s2 = -2*p, s1^2\4 + 2*p, if(!weilok(p,s1,s2), next); my(chi = x^4 - s1*x^3 + s2*x^2 - p*s1*x + p^2); ncls++; my(nJ = subst(chi,x,1), as = -1); for(t=1,#targets, my(T=targets[t], good = (nJ % T[2] == 0)); if(good, for(j=1,#T[3], if(!rank2ok(chi,T[3][j],nJ), good=0; break))); if(good, comp[t]++; if(as<0, as=abssimple(chi)); if(as, simp[t]++; if(#ex[t]<2, ex[t]=concat(ex[t],[[s1,s2,nJ]]))))))); printf("\n########## p = %d : %d Weil-admissible isogeny classes\n", p, ncls); for(t=1,#targets, if(comp[t]==0, badprimes[t]=concat(badprimes[t],[p])); if(simp[t]>0, certprimes[t]=concat(certprimes[t],[p])); printf("  %s  compatible: %5d   ABS-SIMPLE: %5d   %s %s\n", targets[t][1], comp[t], simp[t], if(comp[t]==0, ">>> NO compatible class: forced BAD REDUCTION here", if(simp[t]==0, ">>> compatible but none abs-simple: p can never CERTIFY", "")), if(simp[t]>0 && #ex[t]>0, Str("e.g. ", ex[t]), "")))); printf("\n\n================ SUMMARY over p in %s\n", plist); for(t=1,#targets, printf("  %s  primes with NO compatible class (forced bad reduction): %s\n", targets[t][1], badprimes[t]); printf("  %s  #primes that CAN certify simplicity: %d of %d ; smallest: %s\n", targets[t][1], #certprimes[t], #plist, if(#certprimes[t]>0, certprimes[t][1], "NONE")));}

plist = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293];
census(plist);
quit;
