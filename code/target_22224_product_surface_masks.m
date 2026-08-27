//////////////////////////////////////////////////////////////////////
// New rational surface on the full A(2,2,2,8) cover:
//
//   (a,b,c,d) = (1,r,s,r*s),
//   q0 = 2*t/(1-t^2),
//   s = (1-q0^2*r)/(q0^2-r).
//
// This file computes rigorous projective-r masks on fixed rational t-fibres.
// A residue is retained iff the curve is singular or 3 divides #J(F_p).
//////////////////////////////////////////////////////////////////////
SetColumns(0);
if not assigned PrimeList then
    PrimeList:=[5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
               67,71,73,79,83,89,97,101,103,107,109,113,127,
               131,137,139,149,151,157,163,167,173,179,181,
               191,193,197,199];
elif Type(PrimeList) eq MonStgElt then
    PrimeList:=[StringToInteger(z):z in Split(PrimeList,",")];
end if;
if not assigned output_prefix then
    output_prefix:="results/target_22224_product_surface_masks";
end if;

// <label, numerator(t), denominator(t)>.
fibres:=[<"t1_2",1,2>,<"t2_3",2,3>];
Z:=Integers();

function TupleOnFibre(F,n,d,tn,td)
    t:=F!tn/F!td;
    q:=2*t/(1-t^2);
    A:=q^2;
    // D and N are homogeneous representatives of q^2-r and 1-q^2*r.
    D:=A*d-n;
    N:=d-A*n;
    return [d*D,n*D,d*N,n*N];
end function;

function SmoothTuple(v)
    sq:=[z^2:z in v];
    return &and[z ne 0:z in sq] and #Set(sq) eq 4;
end function;

for p in PrimeList do
    if p in {2,3} then continue; end if;
    F:=GF(p); P<x>:=PolynomialRing(F);
    filename:=Sprintf("%o_p%o.tsv",output_prefix,p);
    out:=Open(filename,"w");
    fprintf out,"fiber\tn\td\tstatus\tjacobian_order\n";
    for fib in fibres do
        label,tn,td:=Explode(fib);
        if F!td eq 0 or (F!td)^2-(F!tn)^2 eq 0 then
            print "PRODUCT_FIBRE_MASK_UNDEFINED","p",p,"fiber",label;
            continue;
        end if;
        boundary:=[]; allowed:=[]; killed:=0;
        pts:=[<z,F!1>:z in F] cat [<F!1,F!0>];
        for pt in pts do
            n,d:=Explode(pt);
            v:=TupleOnFibre(F,n,d,tn,td);
            if not SmoothTuple(v) then
                Append(~boundary,<Z!n,Z!d>);
                fprintf out,"%o\t%o\t%o\tboundary\t0\n",label,Z!n,Z!d;
                continue;
            end if;
            f:=x*&*[x+z^2:z in v];
            nJ:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
            if nJ mod 3 eq 0 then
                Append(~allowed,<Z!n,Z!d>);
                fprintf out,"%o\t%o\t%o\tgood_allowed\t%o\n",
                        label,Z!n,Z!d,nJ;
            else
                killed+:=1;
            end if;
        end for;
        print "PRODUCT_FIBRE_MASK","p",p,"fiber",label,
              "boundary",boundary,"good_allowed",allowed,
              "killed",killed,"total_allowed",#boundary+#allowed;
    end for;
    delete out;
    print "PRODUCT_FIBRE_MASK_FILE",filename;
end for;
