//////////////////////////////////////////////////////////////////////
// Frobenius cycle diagnostics from good finite fibers of the full
// P8-pulled-back independent-3 cover (with L retained).
//
// For every good u in F_p, construct the saturated zero-dimensional
// fiber, take a primitive linear form, and factor its minimal polynomial.
// The factor degrees are Frobenius cycle lengths on that fiber.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned p then p:=7;
elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
assert IsPrime(p) and p notin {2,3,5};

k:=GF(p); PR<q>:=PolynomialRing(k);
tnum:=k!4*(q^2+q-k!6); tden:=q^2+k!6;
enum:=-(k!25/k!3)*tnum^2*tden^2;
eden:=tnum^4-k!25*tnum^2*tden^2+(k!1250/k!3)*tden^4;

procedure AnalyzeFiber(label,ee)
    R<L,U,v>:=PolynomialRing(k,3,"grevlex");
    M:=L^2;
    N:=(k!3*U+k!6*M)/k!2;
    Rc:=(k!3*U^2+k!3*v^2+(k!2/ee-k!15)*M-N^2)/k!2;
    F3:=k!2*v^3+k!2*N*Rc-U^3-k!6*U*v^2-k!22*M;
    F2:=Rc^2+k!2*N*v^3-k!3*U^2*v^2-k!3*v^4
        -(k!1/ee^2-k!15)*M;
    F1:=k!2*Rc*v^3-k!3*U*v^4-(k!2/ee+k!6)*M;
    I:=ideal<R|F3,F2,F1>;
    I:=Saturation(I,ideal<R|L>);
    I:=Saturation(I,ideal<R|v>);
    I:=Saturation(I,ideal<R|U^2-k!4*v^2>);
    if Dimension(I) ne 0 then
        print "FIBER",label,"E",ee,"DIMENSION",Dimension(I);
        return;
    end if;
    A,mp:=quo<R|I>; d:=Dimension(A);
    found:=false;
    for c in [1..6] do
        f:=MinimalPolynomial(mp(L+k!c*U+k!(c*c+1)*v));
        if Degree(f) eq d and GCD(f,Derivative(f)) eq 1 then
            print "FIBER",label,"E",ee,"LENGTH",d,
                  "LINEAR_FORM",c,"CYCLE_DEGREES",
                  Sort([Degree(fe[1]):fe in Factorization(f)]);
            found:=true; break;
        end if;
    end for;
    if not found then print "FIBER",label,"E",ee,"LENGTH",d,
                            "NO_SQUAREFREE_PRIMITIVE_FORM"; end if;
end procedure;

print "CONTACT6_M612_P8_CORE_FIBER_CYCLES_ROOT";
print "PRIME",p;
for q0 in k do
    de:=Evaluate(eden,q0);
    if de eq 0 then print "PARAMETER_POLE",Integers()!q0; continue; end if;
    ee:=Evaluate(enum,q0)/de;
    if ee eq 0 or k!1+k!2*ee eq 0 then
        print "PARAMETER_BOUNDARY",Integers()!q0,"E",ee; continue;
    end if;
    AnalyzeFiber(Integers()!q0,ee);
end for;
einf:=LeadingCoefficient(enum)/LeadingCoefficient(eden);
if einf ne 0 and k!1+k!2*einf ne 0 then AnalyzeFiber("infinity",einf); end if;
print "CONTACT6_M612_P8_CORE_FIBER_CYCLES_ROOT_DONE";
quit;
