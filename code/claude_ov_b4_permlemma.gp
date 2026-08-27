/* claude_ov_b4_permlemma.gp -- Lane 4 (route B4): the PERMUTATION LEMMA that
   closes the transverse (D6) component of the Richelot-source locus.

   Lemma.  Let f = c*G1*G2*G3 with Gi monic quadratics whose SET is Galois
   stable, Delta = det[coeffs], H_i = G_j' G_k - G_j G_k' ((i,j,k) cyclic),
   codomain f' = -H1 H2 H3 / Delta.  Then for every permutation sigma of the
   indices, H_i |-> eps_i H_{sigma(i)} with eps_i = +-1 and Delta |-> sgn(sigma)
   Delta, so Gal permutes the THREE QUADRATIC FACTORS of f' by exactly the same
   permutation representation rho: Gal -> S3 that it uses on {G1,G2,G3}.

   Consequence.  If rho is TRANSITIVE, every Galois orbit on the 6 roots of f'
   has size divisible by 3, so the rational factor type of f' is [6] or [3,3],
   and in both cases the number of even subsets is 2, i.e. the codomain has
   2-RANK 0.  If rho is NOT transitive (image trivial or of order 2) then some
   G_i is Galois stable, i.e. f has a RATIONAL quadratic factor and the member
   lies on the quadratic-factor incidence variety.

   This script verifies the lemma symbolically with generic G1,G2,G3.
   Usage: gp -q code/claude_ov_b4_permlemma.gp
*/
x='x;
G(a,b) = x^2 + a*x + b;
mkdel(g1,g2,g3) = matdet([polcoef(g1,2,x),polcoef(g1,1,x),polcoef(g1,0,x); polcoef(g2,2,x),polcoef(g2,1,x),polcoef(g2,0,x); polcoef(g3,2,x),polcoef(g3,1,x),polcoef(g3,0,x)]);
H(gj,gk) = deriv(gj,x)*gk - gj*deriv(gk,x);
Hs(g1,g2,g3) = [H(g2,g3), H(g3,g1), H(g1,g2)];
{
  my(a1='a1,b1='b1,a2='a2,b2='b2,a3='a3,b3='b3, g1,g2,g3, HH, DD, ok);
  g1 = G(a1,b1); g2 = G(a2,b2); g3 = G(a3,b3);
  HH = Hs(g1,g2,g3); DD = mkdel(g1,g2,g3);
  /* transposition (2 3) */
  my(HHt = Hs(g1,g3,g2), DDt = mkdel(g1,g3,g2));
  print("swap(2,3): H1 -> -H1 ? ", HHt[1] == -HH[1]);
  print("swap(2,3): H2 -> -H3 ? ", HHt[2] == -HH[3]);
  print("swap(2,3): H3 -> -H2 ? ", HHt[3] == -HH[2]);
  print("swap(2,3): Delta -> -Delta ? ", DDt == -DD);
  print("swap(2,3): H1H2H3/Delta INVARIANT ? ",
        HHt[1]*HHt[2]*HHt[3]/DDt == HH[1]*HH[2]*HH[3]/DD);
  /* 3-cycle (1 2 3): (g1,g2,g3) -> (g3,g1,g2) */
  my(HHc = Hs(g3,g1,g2), DDc = mkdel(g3,g1,g2));
  print("cycle: {H} permuted, product/Delta INVARIANT ? ",
        HHc[1]*HHc[2]*HHc[3]/DDc == HH[1]*HH[2]*HH[3]/DD);
  /* Hs(g3,g1,g2) returns [H(g1,g2), H(g2,g3), H(g3,g1)] = [H3, H1, H2] */
  print("cycle: the H-triple is the cyclic shift of the old one ? ",
        HHc == [HH[3], HH[1], HH[2]]);
  print("cycle: Delta invariant ? ", DDc == DD);
  /* each H_i is a quadratic in x */
  print("deg H_i = ", [poldegree(h,x) | h <- HH]);
}
/* even-subset count for the two possible transitive factor types */
tworank(degs) = { my(k=#degs, ne=0); for(m=0,2^k-1, my(s=0); for(i=1,k, if(bittest(m,i-1), s+=degs[i])); if(s%2==0, ne++)); round(log(ne)/log(2))-1; };
print("2-rank of factor type [6]   = ", tworank([6]));
print("2-rank of factor type [3,3] = ", tworank([3,3]));
print("PERMLEMMA_DONE");
quit;
