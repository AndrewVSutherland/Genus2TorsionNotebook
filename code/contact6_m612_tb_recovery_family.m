//////////////////////////////////////////////////////////////////////
// The rational one-parameter recovery-boundary family on the compact
// T_B-halving cover.
//
// Modes:
//   summary  exact two-branch formulas and specialized cubic-contact cover;
//   masks    exact finite extra-3 masks on the u-line.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode:="summary"; end if;
if not assigned prime_bound then prime_bound:=97;
elif Type(prime_bound) eq MonStgElt then
    prime_bound:=StringToInteger(prime_bound);
end if;

Q:=Rationals(); Z:=Integers();

function Has212(invs)
    return #[n:n in invs|(Z!n) mod 2 eq 0] ge 2
       and #[n:n in invs|(Z!n) mod 12 eq 0] ge 1;
end function;

function Has612(invs)
    return #[n:n in invs|(Z!n) mod 6 eq 0] ge 2
       and #[n:n in invs|(Z!n) mod 12 eq 0] ge 1;
end function;

function Primitive(poly)
    if poly eq 0 then return poly; end if;
    den:=LCM([Denominator(c):c in Coefficients(poly)]);
    vals:=[Z!(den*c):c in Coefficients(poly)];
    cont:=GCD(vals); if cont eq 0 then cont:=1; end if;
    return Parent(poly)!((Q!den/Q!Abs(cont))*poly);
end function;

if mode eq "summary" then
    R<u,M,U,v>:=PolynomialRing(Q,4,"grevlex");
    Kf:=FieldOfFractions(R);
    uu:=Kf!u;
    D:=(uu^2-1)*(uu^2-9);
    s:=-(3*uu^4-14*uu^2+27)/(2*D);
    y:=(uu^4-9)/D;
    w:=-8*uu*(uu^2-3)/D;
    r:=s*(2*s+3);
    m:=(s+1)*w;
    KK:=m^2;
    b:=2*s^2-3;
    abase:=16*s^4+36*s^3+18*s^2-4*s-3;

    print "CONTACT6_M612_TB_RECOVERY_FAMILY";
    print "D",D;
    print "s",s;
    print "y",y;
    print "w",w;
    print "IDENTITIES",
          y^2 eq (s+1)*(s-1/2),
          w^2 eq (2*s-1)*(2*s+3),
          KK eq (s+1)^2*(2*s-1)*(2*s+3);
    print "r",r;
    print "m",m;
    print "K",KK;
    print "b",b;

    covers:=[];
    for eps in [-1,1] do
        a:=abase+eps*16*s*(s+1)^2*y;
        A3:=a-3+4*s^2*r-2*KK;
        H1:=s*((a-3)*r^2+4*r-KK*(a+3))-r*A3;
        H2:=8*s^2*(2*s^2*r^2+2*(a-3)*r+2
                    -KK*(2*s^2-6))-A3^2-32*s^3*r;
        assert H1 eq 0 and H2 eq 0;
        print "BRANCH",eps,"a",a;

        c1:=2*a+6;
        c2:=a^2+2*b-15;
        c3:=2*a*b+22;
        c4:=2*a+b^2-15;
        c5:=2*b+6;
        BB:=c5*M+3*U;
        Delta:=4*c4*M+12*(U^2+v^2)-BB^2;
        F3:=BB*Delta+16*v^3-8*c3*M-8*U^3-48*U*v^2;
        F2:=Delta^2+64*BB*v^3-64*c2*M
             -192*(U^2*v^2+v^4);
        F1:=Delta*v^3-4*c1*M-12*U*v^4;
        cores:=[Primitive(R!Numerator(F)):F in [F1,F2,F3]];
        print " CUBIC_CONTACT_SHAPES",
              [<TotalDegree(F),Degree(F,u),Degree(F,M),Degree(F,U),
                Degree(F,v),#Terms(F)>:F in cores];
        print " CUBIC_CONTACT_GCD",GCD(GCD(cores[1],cores[2]),cores[3]);
        Append(~covers,cores);
    end for;
    print "NONBOUNDARY",
          "u*(u^2-1)*(u^2-9)*M*v*(v^3-1)*(U^2-4v^2)*Disc(f) != 0";
    print "CONTACT6_M612_TB_RECOVERY_FAMILY_DONE";
    quit;
end if;

if mode ne "masks" then error "mode must be summary or masks"; end if;

machine:=0;
if assigned output_machine then
    machine:=Open(output_machine,"w");
end if;
print "CONTACT6_M612_TB_RECOVERY_MASKS","prime_bound",prime_bound;
for p in [ell:ell in PrimesUpTo(prime_bound)|ell notin {2,3}] do
    k:=GF(p); Pk<X>:=PolynomialRing(k);
    for eps in [-1,1] do
        allowed:=[]; bad:=[]; base212:=0; good:=0;
        for u in k do
            D:=(u^2-1)*(u^2-9);
            if D eq 0 then Append(~bad,Z!u); continue; end if;
            s:=-(3*u^4-14*u^2+27)/(2*D);
            y:=(u^4-9)/D;
            w:=-8*u*(u^2-3)/D;
            r:=s*(2*s+3);
            m:=(s+1)*w;
            if s*r*m eq 0 then Append(~bad,Z!u); continue; end if;
            b:=2*s^2-3;
            abase:=16*s^4+36*s^3+18*s^2-4*s-3;
            a:=abase+eps*16*s*(s+1)^2*y;
            h:=1+a*X+b*X^2+X^3;
            f:=h^2-(X-1)^6;
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                Append(~bad,Z!u); continue;
            end if;
            good+:=1;
            AG,phi:=AbelianGroup(Jacobian(HyperellipticCurve(f)));
            inv:=Invariants(AG);
            if Has212(inv) then base212+:=1;
            else
                print "BASE212_ASSERT_FAIL",p,eps,Z!u,inv;
            end if;
            if Has612(inv) then Append(~allowed,Z!u); end if;
        end for;
        print "MASK","p",p,"branch",eps,"good",good,
              "base212",base212,"allowed",#allowed,"bad",#bad;
        print " ALLOWED",allowed;
        print " BAD",bad;
        if assigned output_machine then
            pass:=Sort(Setseq(Seqset(allowed cat bad)));
            fprintf machine,"%o %o %o",p,eps,#pass;
            for z in pass do fprintf machine," %o",z; end for;
            fprintf machine,"\n";
        end if;
    end for;
end for;
if assigned output_machine then delete machine; end if;
print "CONTACT6_M612_TB_RECOVERY_MASKS_DONE";
quit;
