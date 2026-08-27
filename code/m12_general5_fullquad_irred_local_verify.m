//////////////////////////////////////////////////////////////////////
// Exact finite-field regression for the nonsplit quadratic-B chart.
//
// These are the two open nonsplit-B points returned at p=7 by
// m12_general5_fullquad_irred.py.  Besides checking the norm identity,
// this script constructs the Mumford 5-class and the marked M(12)
// 12-class and verifies that their sum has exact order 60.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
k := GF(7);
P<x> := PolynomialRing(k);

rows := [
    // <b,w,r,n,s,t,a,tau,lambda,q,A>
    <k!3,k!5,k!1,k!4,k!4,k!2,k!3,k!4,k!2,
     x^2+x,
     x^5+5*x^4+5*x^3+6*x^2+3*x+4>,
    <k!5,k!6,k!6,k!3,k!2,k!3,k!3,k!2,k!3,
     x^2+x,
     x^5+x^4+4*x^3+x^2+5>
];

for row in rows do
    b,w,r,n,s,t,a,tau,lambda,q,A := Explode(row);
    B := x^2+r*x+n;
    assert not IsSquare(r^2-4*n);
    assert IsSquarefree(B) and IsSquarefree(q) and GCD(B,q) eq 1;

    L := b+(2*b-1)*x;
    H := x+w*(1+b*x);
    F := P!(L*(L*H^2+4*b*(1+x)^2*(w*L-x^2)));
    assert Degree(F) eq 5 and IsSquarefree(F);
    assert lambda^2 eq tau;
    assert A^2-(lambda*B)^2*F eq q^5;
    assert GCD(q,F) eq 1;

    v := (-A*InverseMod(lambda*B,q)) mod q;
    assert (v^2-F) mod q eq 0;
    C := HyperellipticCurve(F);
    J := Jacobian(C);
    D5 := J![q,v];

    eta := (1-b)*(w*(1-b)-1);
    assert Evaluate(F,-1) eq eta^2;
    D12 := J![x+1,eta];

    assert Order(D5) eq 5;
    assert Order(D12) eq 12;
    assert Order(D5+D12) eq 60;
    print "IRRED_LOCAL_CONTROL", "b",b,"w",w,"B",B,
          "q",q,"v",v,"#J",#J,
          "orders",<Order(D5),Order(D12),Order(D5+D12)>;
end for;

print "M12_GENERAL5_FULLQUAD_IRRED_LOCAL_VERIFY_PASS";
quit;
