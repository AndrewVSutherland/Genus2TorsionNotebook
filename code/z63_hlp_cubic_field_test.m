//////////////////////////////////////////////////////////////////////
// Sufficient cubic-field test for the low-height HLP 7-by-9
// squareclass candidates.  Equal discriminant squareclass is only a
// necessary filter; isomorphic 2-torsion cubic fields are required.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemMB then MemMB:=180; end if;
if Type(MemMB) eq MonStgElt then MemMB:=StringToInteger(MemMB); end if;
SetMemoryLimit(MemMB*10^6);

Q:=Rationals(); P<x>:=PolynomialRing(Q);

function TwoDivisionCubic(N,t)
    if N eq 7 then
        b:=t^3-t^2; c:=t^2-t;
    elif N eq 9 then
        b:=t^5-2*t^4+2*t^3-t^2; c:=t^3-t^2;
    else
        error "Only N=7,9 are supported";
    end if;
    return 4*x^3-4*b*x^2+((1-c)*x-b)^2;
end function;

function TateCurve(N,t)
    if N eq 7 then
        b:=t^3-t^2; c:=t^2-t;
    elif N eq 9 then
        b:=t^5-2*t^4+2*t^3-t^2; c:=t^3-t^2;
    end if;
    return EllipticCurve([1-c,-b,-b,Q!0,Q!0]);
end function;

T7:=[Q|-16/3,19/16,3/19,-2,3/2,1/3];
T9:=[Q|4,-1/3,3/4,-1,2,1/2];

print "Z63_HLP_CUBIC_FIELD_TEST";
J7:=[jInvariant(TateCurve(7,t)):t in T7];
J9:=[jInvariant(TateCurve(9,u)):u in T9];
assert #SequenceToSet(J7[1..3]) eq 1 and #SequenceToSet(J7[4..6]) eq 1;
assert #SequenceToSet(J9[1..3]) eq 1 and #SequenceToSet(J9[4..6]) eq 1;
// Certify that each three-element block is one Q-isomorphism orbit, not
// merely a set of twists with equal j-invariant.
assert &and[IsIsomorphic(TateCurve(7,T7[1]),TateCurve(7,t))
            : t in T7[2..3]];
assert &and[IsIsomorphic(TateCurve(7,T7[4]),TateCurve(7,t))
            : t in T7[5..6]];
assert &and[IsIsomorphic(TateCurve(9,T9[1]),TateCurve(9,u))
            : u in T9[2..3]];
assert &and[IsIsomorphic(TateCurve(9,T9[4]),TateCurve(9,u))
            : u in T9[5..6]];
print "J_ORBITS_7",J7[1],J7[4];
print "J_ORBITS_9",J9[1],J9[4];
hits:=[];
for t in T7 do
    g7:=TwoDivisionCubic(7,t);
    assert IsIrreducible(g7);
    K7<a>:=NumberField(g7);
    for u in T9 do
        // The first three values in each list have squareclass -741;
        // the last three have squareclass -6.  Skip cross-class pairs.
        if (Index(T7,t) le 3) ne (Index(T9,u) le 3) then continue; end if;
        g9:=TwoDivisionCubic(9,u);
        assert IsIrreducible(g9);
        K9<b>:=NumberField(g9);
        iso:=IsIsomorphic(K7,K9);
        print "PAIR",t,u,"field_disc",Discriminant(MaximalOrder(K7)),
              Discriminant(MaximalOrder(K9)),"isomorphic",iso;
        if iso then Append(~hits,<t,u>); end if;
    end for;
end for;
print "ISOMORPHIC_PAIRS",#hits,hits;
assert #hits eq 9;
assert <-16/3,Q!4> in hits;
print "Z63_HLP_CUBIC_FIELD_TEST_DONE";
quit;
