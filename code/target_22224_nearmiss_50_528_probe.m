// Exact and twist-local probe of the first transverse all-contact-mask near miss.
SetColumns(0);
SetLogFile("results/target_22224_nearmiss_50_528_probe.log" : Overwrite := true);
Q:=Rationals(); Z:=Integers(); P<x>:=PolynomialRing(Q);
v:=[50,528,-726,-891];
f:=x*&*[x+(Q!z)^2:z in v];
C:=HyperellipticCurve(f); J:=Jacobian(C);
T:=TorsionSubgroup(J);
print "NEARMISS_START","tuple",v,"f",f,"torsion",Invariants(T);
print "LOCAL_RESIDUES_11",[z mod 11:z in v],"squares",[z^2 mod 11:z in v];
print "LOCAL_RESIDUES_13",[z mod 13:z in v],"squares",[z^2 mod 13:z in v];
for p in [17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101] do
    k:=GF(p); R<X>:=PolynomialRing(k);
    fp:=X*&*[X+(k!(z mod p))^2:z in v];
    if Discriminant(fp) eq 0 then continue; end if;
    eta:=k!2; while IsSquare(eta) do eta+:=1; end while;
    np:=#Jacobian(HyperellipticCurve(fp));
    nm:=#Jacobian(HyperellipticCurve(eta*fp));
    print "TWIST_LOCAL",p,np,nm,"div192",np mod 192 eq 0,nm mod 192 eq 0;
    if np mod 192 ne 0 and nm mod 192 ne 0 then
        print "TWIST_FAMILY_OBSTRUCTION",p,np,nm;
        break;
    end if;
end for;
print "NEARMISS_DONE";
quit;
