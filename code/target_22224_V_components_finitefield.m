// Independent finite-field check for the V component analysis.
// At cover-allowed finite residues we report smoothness and #J(F_p).
SetColumns(0);
SetLogFile("results/target_22224_V_components_finitefield.log":Overwrite:=true);
procedure Check(name,p,vals,Ts)
    F:=GF(p); P<x>:=PolynomialRing(F);
    for t in Ts do
        base:=[F!vals[1],F!vals[2],F!vals[3],F!vals[4]*(F!t)^2];
        f:=x*&*[x+a^2:a in base];
        smooth:=Degree(f) eq 5 and Discriminant(f) ne 0;
        if smooth then
            nJ:=Integers()!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
            print "FINITE_FIELD",name,"p",p,"T",t,"base",base,
                  "smooth",smooth,"Jorder",nJ,"mod3",nJ mod 3;
        else
            print "FINITE_FIELD",name,"p",p,"T",t,"base",base,
                  "smooth",smooth;
        end if;
    end for;
end procedure;
Check("V1",13,[-2178,2420,9075,-1470],[1,5,8,12]);
Check("V2",11,[-1458,2268,7938,-2023],[1,5,6,10]);
Check("V2",13,[-1458,2268,7938,-2023],[1,4,9,12]);
UnsetLogFile();
quit;
