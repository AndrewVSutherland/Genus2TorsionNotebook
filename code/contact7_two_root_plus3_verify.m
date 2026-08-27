//////////////////////////////////////////////////////////////////////
// Independent Magma verification for the two-root contact-7 plus-3 sieve.
//
// Signed parameters absorb the original signs:
//     r=1-u^2, w=1-v^2, h(r)=u^7, h(w)=v^7.
//
// Modes:
//   magma -b mode:=finite prime:=11 code/contact7_two_root_plus3_verify.m
//   magma -b mode:=sample u:=2 v:=3 code/contact7_two_root_plus3_verify.m
//   magma -b mode:=verify \
//       input_file:="data/contact7_two_root_plus3_h100_survivors.m" \
//       code/contact7_two_root_plus3_verify.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "sample"; end if;
if not assigned prime then prime := 11;
elif Type(prime) eq MonStgElt then prime := StringToInteger(prime); end if;
if not assigned u then u := 2; elif Type(u) eq MonStgElt then u := StringToRational(u); end if;
if not assigned v then v := 3; elif Type(v) eq MonStgElt then v := StringToRational(v); end if;
if not assigned input_file then
    input_file := "data/contact7_two_root_plus3_h100_survivors.m";
end if;
if not assigned max_exact then max_exact := 10000;
elif Type(max_exact) eq MonStgElt then max_exact := StringToInteger(max_exact); end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function RequiredPrimeToP(n,p)
    while n mod p eq 0 do n div:= p; end while;
    return n;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L,Denominator(Coefficient(f,i)));
    end for;
    return P!(L^2*f),L;
end function;

function TwoRootModel(uu,vv)
    uu:=Q!uu; vv:=Q!vv;
    r:=1-uu^2; w:=1-vv^2;
    if uu eq 0 or vv eq 0 then
        return false,"singular contact collision",P!0,P!0,Q!0,Q!0,r,w;
    end if;
    if uu eq -1 or vv eq -1 or uu eq -vv then
        return false,"incompatible parameter pole",P!0,P!0,Q!0,Q!0,r,w;
    end if;
    if uu eq vv then
        return false,"coincident-root singularity",P!0,P!0,Q!0,Q!0,r,w;
    end if;
    Auu:=(2*uu^5+4*uu^4+6*uu^3+8*uu^2+10*uu+5)/(2*(uu+1)^2);
    Avv:=(2*vv^5+4*vv^4+6*vv^3+8*vv^2+10*vv+5)/(2*(vv+1)^2);
    b:=(Avv-Auu)/(uu^2-vv^2);
    a:=Auu-b*r;
    h:=1-(Q!7/2)*x+a*x^2+b*x^3;
    f:=ExactQuotient(h^2+(x-1)^7,x^2);
    assert Evaluate(h,r) eq uu^7 and Evaluate(h,w) eq vv^7;
    assert Evaluate(f,r) eq 0 and Evaluate(f,w) eq 0;
    assert h^2-x^2*f eq -(x-1)^7;
    if Evaluate(h,1) eq 0 then
        return false,"contact boundary",f,h,a,b,r,w;
    end if;
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false,"singular",f,h,a,b,r,w;
    end if;
    return true,"",f,h,a,b,r,w;
end function;

function TwoRootModelFinite(F,uu,vv)
    PF<X>:=PolynomialRing(F);
    r:=1-uu^2; w:=1-vv^2;
    if uu eq 0 or vv eq 0 or uu eq -F!1 or vv eq -F!1 or
       uu eq vv or uu eq -vv then
        return false,"chart",PF!0,PF!0,F!0,F!0,r,w;
    end if;
    Auu:=(2*uu^5+4*uu^4+6*uu^3+8*uu^2+10*uu+5)/(2*(uu+1)^2);
    Avv:=(2*vv^5+4*vv^4+6*vv^3+8*vv^2+10*vv+5)/(2*(vv+1)^2);
    b:=(Avv-Auu)/(uu^2-vv^2);
    a:=Auu-b*r;
    h:=1-(F!7/F!2)*X+a*X^2+b*X^3;
    if Evaluate(h,F!1) eq 0 then
        return false,"contact",PF!0,h,a,b,r,w;
    end if;
    f:=ExactQuotient(h^2+(X-1)^7,X^2);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false,"bad",f,h,a,b,r,w;
    end if;
    return true,"",f,h,a,b,r,w;
end function;

function FrobeniusPolynomial(C,p)
    ef:=EulerFactor(C,p);
    d:=Degree(ef);
    return &+[Q!Coefficient(ef,i)*x^(d-i):i in [0..d]];
end function;

function SimplicityCertificate(C)
    qfound:=false; qp:=0; qPhi:=P!0; qord:=0; qdesc:="";
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149] do
        try
            Phi:=FrobeniusPolynomial(C,p);
            fac:=Factorization(Phi);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                // Irreducibility already certifies Q-simplicity.  The D4
                // test is also recorded when it supplies the standard
                // geometric-simplicity certificate for an abelian surface.
                geometric:=false; desc:=""; ord:=0;
                try
                    G:=GaloisGroup(Phi);
                    ord:=Order(G);
                    try desc:=TransitiveGroupDescription(G); catch e2 desc:=""; end try;
                    geometric := ord eq 8 and desc eq "D(4)";
                catch e3
                    geometric:=false;
                end try;
                if not qfound then
                    qfound:=true; qp:=p; qPhi:=Phi; qord:=ord; qdesc:=desc;
                end if;
                if geometric then return true,true,p,Phi,ord,desc; end if;
            end if;
        catch e
            continue;
        end try;
    end for;
    if qfound then return true,false,qp,qPhi,qord,qdesc; end if;
    return false,false,0,P!0,0,"";
end function;

function PassesAdditionalReductions(f)
    // These begin beyond the C++ masks.  Reduction is used only when the
    // displayed rational polynomial is p-integral, degree 5, and squarefree;
    // all denominator and bad-reduction primes are skipped conservatively.
    plist:=[89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,
            173,179,181,191,193,197,199,211,223,227,229,233,239,241,251];
    used:=[];
    for p in plist do
        try
            fp:=ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            n:=Integers()!#Jacobian(HyperellipticCurve(fp));
            Append(~used,<p,n>);
            req:=RequiredPrimeToP(84,p);
            if n mod req ne 0 then return false,p,n,used; end if;
        catch e
            continue;
        end try;
    end for;
    return true,0,0,used;
end function;

procedure VerifyOne(uu,vv,DoTorsion)
    ok,msg,f,h,a,b,r,w:=TwoRootModel(uu,vv);
    print "PARAMETERS",uu,vv,"ok",ok,"message",msg;
    if not ok then return; end if;
    print "a",a,"b",b,"r",r,"w",w;
    print "f",f;
    print "factorization",Factorization(f);
    C:=HyperellipticCurve(f); J:=Jacobian(C);
    D7:=J![x-1,Evaluate(h,1)]; Dr:=J![x-r,0]; Dw:=J![x-w,0];
    print "visible_orders",Order(D7),Order(Dr),Order(Dw),Order(Dr+Dw);
    assert Order(D7) eq 7;
    assert Order(Dr) eq 2 and Order(Dw) eq 2 and Order(Dr+Dw) eq 2 and Dr ne Dw;
    if not DoTorsion then return; end if;
    fI,L:=IntegralModel(f);
    CI:=HyperellipticCurve(fI); JI:=Jacobian(CI);
    // Recheck the marked classes after y -> L*y.
    D7I:=JI![x-1,L*Evaluate(h,1)];
    DrI:=JI![x-r,0]; DwI:=JI![x-w,0];
    assert Order(D7I) eq 7;
    assert Order(DrI) eq 2 and Order(DwI) eq 2 and Order(DrI+DwI) eq 2;
    G,phi:=TorsionSubgroup(JI);
    invs:=Invariants(G);
    qsimp,gsimp,p,Phi,gord,gdesc:=SimplicityCertificate(CI);
    print "EXACT_TORSION",invs,"order",#G;
    print "SIMPLICITY", "Q_simple",qsimp,"geometric_D4",gsimp,
          "prime",p,"Galois_order",gord,"description",gdesc;
    if qsimp then print "Frobenius",Phi; end if;
    if #G gt 80 then
        print "RECORD_HIT",uu,vv,"torsion",invs,"order",#G;
        print "integral_model",fI;
    end if;
end procedure;

if mode eq "finite" then
    p:=prime; F:=GF(p);
    total:=0; chart:=0; good:=0; pass:=0; fail:=0; base_anomaly:=0;
    pass_samples:=[]; fail_samples:=[];
    req:=RequiredPrimeToP(84,p); base:=RequiredPrimeToP(28,p);
    for uu in F do for vv in F do
        total+:=1;
        ok,msg,f,h,a,b,r,w:=TwoRootModelFinite(F,uu,vv);
        if msg ne "chart" then chart+:=1; end if;
        if not ok then continue; end if;
        good+:=1;
        n:=Integers()!#Jacobian(HyperellipticCurve(f));
        if n mod base ne 0 then base_anomaly+:=1; end if;
        rec:=<Integers()!uu,Integers()!vv,Integers()!a,Integers()!b,n>;
        if n mod req eq 0 then
            pass+:=1; if #pass_samples lt 5 then Append(~pass_samples,rec); end if;
        else
            fail+:=1; if #fail_samples lt 5 then Append(~fail_samples,rec); end if;
        end if;
    end for; end for;
    print "FINITE","prime",p,"total",total,"chart",chart,"good",good,
          "pass",pass,"fail",fail,"unknown",total-good,
          "base_anomaly",base_anomaly,"required",req;
    print "PASS_SAMPLES",pass_samples;
    print "FAIL_SAMPLES",fail_samples;
elif mode eq "sample" then
    VerifyOne(Q!u,Q!v,true);
elif mode eq "verify" then
    candidates:=eval Read(input_file);
    print "VERIFY_FILE",input_file,"candidates",#candidates,"max_exact",max_exact;
    done:=0; local_survivors:=0; first_kill:=AssociativeArray();
    for rec in candidates do
        if done ge max_exact then break; end if;
        uu:=Q!rec[1]/rec[2]; vv:=Q!rec[3]/rec[4];
        done+:=1;
        ok,msg,f,h,a,b,r,w:=TwoRootModel(uu,vv);
        if not ok then
            print "INVALID_CANDIDATE",uu,vv,msg;
            continue;
        end if;
        pass,pbad,nbad,used:=PassesAdditionalReductions(f);
        if not pass then
            if IsDefined(first_kill,pbad) then first_kill[pbad]+:=1;
            else first_kill[pbad]:=1; end if;
            print "LOCAL_KILL",uu,vv,"prime",pbad,"Jacobian_order",nbad;
            continue;
        end if;
        local_survivors+:=1;
        print "ADDITIONAL_LOCAL_SURVIVOR",uu,vv,"used",used;
        VerifyOne(uu,vv,true);
    end for;
    print "VERIFY_DONE",done,"additional_local_survivors",local_survivors;
    print "FIRST_KILL";
    for p in Sort([k:k in Keys(first_kill)]) do print p,first_kill[p]; end for;
else
    print "unknown mode",mode;
end if;

quit;
