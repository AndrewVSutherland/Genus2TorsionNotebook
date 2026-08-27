//////////////////////////////////////////////////////////////////////
// Strict first blow-ups of the forced mod-5 contact-core boundary.
//
// On each affine base support g=0, every F_5 point of the cubic-contact
// core with M=L^2 lies on one of
//
//   L=0, v=0, U-2v=0, U+2v=0.
//
// For every pair (g,D) this script checks the two ordinary blow-up charts
//
//   chart 1: e=q*z, D=q,       direction [e:D]=[z:1],
//   chart 2: e=q,   D=q*z,     direction [e:D]=[1:z].
//
// IMPORTANT: the proper transform is Saturation(<F1,F2,F3>,<q>), not
// the ideal obtained by stripping a q power from each chosen generator.
// Global saturation is very expensive on these charts.  This script uses
// the following exact local certificate instead.  If the 3 x 4 Jacobian in
// (t,z,X,Y) has rank 3 at q=0, the raw pullback is smooth there and the
// implicit-function theorem gives a component on which q is free.  It is
// therefore not contained in q=0, saturation changes nothing locally, and
// the point is a genuine Hensel-smooth point of the strict transform.
// Rank-deficient raw points are printed as UNRESOLVED and are not claimed.
//
// Each deep chart is built in its own five-variable ring.  Thus there is
// no unused L, U, or v coordinate multiplying point counts by 5.  Points
// with z != 0 are the balanced cones (both e and D have the same leading
// order); the reciprocal overlap of the two charts is deduplicated.
//
// A cone is called STRICT_SMOOTH_CERTIFIED only under this exact rank-3
// local criterion.  Merely solving a first-order inhomogeneous equation at
// rank < 3 is not enough.
//
// For the DB and DC base supports the relevant dual square cover is
//
//   DB=-8e=W^2, respectively DC=-8e=W^2.
//
// On a balanced chart this eliminates q as a unit times W^2.  At W=0
// the pulled-back Jacobian is therefore the core Jacobian in (t,z,X,Y),
// with its q column removed.  DUAL_SQUARE_SMOOTH records exact rank 3 of
// that matrix.  This first weighted layer incorporates both the even
// valuation and the constrained leading unit.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
k:=GF(5);

function CoreEquations(a,b,l,u,w)
    M:=l^2;
    c1:=2*a+6; c2:=a^2+2*b-15; c3:=2*a*b+22;
    c4:=2*a+b^2-15; c5:=2*b+6;
    B3:=c5*M+3*u;
    Delta3:=4*c4*M+12*(u^2+w^2)-B3^2;
    F3:=B3*Delta3+16*w^3-8*c3*M-8*u^3-48*u*w^2;
    F2:=Delta3^2+64*B3*w^3-64*c2*M
        -192*(u^2*w^2+w^4);
    F1:=Delta3*w^3-4*c1*M-12*u*w^4;
    return [F1,F2,F3];
end function;

function BaseParameters(label,e,t)
    inv8:=(k!8)^-1;
    if label eq "b3" then
        return t,-3+e;
    elif label eq "a3" then
        return -3+e,t;
    elif label eq "rr" then
        return t,-t-2+e;
    elif label eq "DB" then
        return t,inv8*(t-3)^2-3+e;
    elif label eq "DC" then
        return inv8*(t-3)^2-3+e,t;
    end if;
    error "unknown base chart",label;
end function;

// R is always k[q,t,z,X,Y].  X,Y are the two surviving core coordinates.
// The deep coordinate itself is q or q*z and is not retained as a dummy
// residue coordinate on the exceptional fiber.
function BlowupExpressions(deep,chart,q,z,X,Y)
    e:=chart eq 1 select q*z else q;
    d:=chart eq 1 select q else q*z;
    if deep eq "L" then
        return e,d,X,Y;             // (e,L,U,v), X=U, Y=v
    elif deep eq "v" then
        return e,X,Y,d;             // (e,L,U,v), X=L, Y=U
    elif deep eq "Up" then
        return e,X,2*Y+d,Y;         // X=L, Y=v, D=U-2v
    elif deep eq "Um" then
        return e,X,-2*Y+d,Y;        // X=L, Y=v, D=U+2v
    end if;
    error "unknown deep divisor",deep;
end function;

function ConeOpen(deep,x0,y0)
    if deep eq "L" then
        return y0 ne 0 and x0^2-4*y0^2 ne 0;
    elif deep eq "v" then
        return x0 ne 0 and y0 ne 0;
    else
        return x0 ne 0 and y0 ne 0;
    end if;
end function;

function ConeOpenPolynomial(deep,X,Y)
    if deep eq "L" then return Y*(X^2-4*Y^2); end if;
    return X*Y;
end function;

function CountRank(ranks,r)
    return #[x:x in ranks|x eq r];
end function;

// Canonical projective direction.  The chart-1 direction is [z:1].
// The chart-2 direction [1:z] is [1/z:1] on the balanced overlap.
function Direction(chart,z0)
    if chart eq 1 then return z0,k!1; end if;
    if z0 eq 0 then return k!1,k!0; end if;
    return z0^-1,k!1;
end function;

procedure AddCone(~keys,~ranks,~square_ranks,key,r,sr)
    for i in [1..#keys] do
        if keys[i] eq key then
            assert ranks[i] eq r;
            assert square_ranks[i] eq sr;
            return;
        end if;
    end for;
    Append(~keys,key);
    Append(~ranks,r);
    Append(~square_ranks,sr);
end procedure;

labels:=["b3","a3","rr","DB","DC"];
deeps:=["L","v","Up","Um"];

print "CONTACT6_M612_AFFINE_FIRST_BLOWUPS_MOD5_STRICT";
print "FIELDS: pair, raw balanced cones, rank distribution, certified strict-smooth,";
print "        dual_square_smooth_rank3 (DB/DC only), smooth_cone_keys";
print "CONE_KEY=<t,[e:D]_0,[e:D]_1,X,Y>; X,Y meanings depend on DEEP";
print "  L: X=U,Y=v; v: X=L,Y=U; Up/Um: X=L,Y=v";

grand_cones:=0;
grand_smooth:=0;
grand_square_smooth:=0;

for label in labels do
    for deep in deeps do
        cone_keys:=[];
        cone_ranks:=[];
        cone_square_ranks:=[];

        for chart in [1,2] do
            R<q,t,z,X,Y>:=PolynomialRing(k,5,"grevlex");
            e,l,u,w:=BlowupExpressions(deep,chart,q,z,X,Y);
            a,b:=BaseParameters(label,e,t);
            raw:=CoreEquations(a,b,l,u,w);

            eqs:=raw;
            derivs:=[[Derivative(f,j):j in [1..5]]:f in eqs];
            balanced:=0;
            full_ranks:=[];
            fiber_ranks:=[];
            square_balanced:=0;
            square_balanced_smooth:=0;

            for t0 in k do for z0 in [z0:z0 in k|z0 ne 0] do
                for x0 in k do for y0 in k do
                if not ConeOpen(deep,x0,y0) then continue; end if;
                pt:=[k!0,t0,z0,x0,y0];
                if not &and[Evaluate(f,pt) eq 0:f in eqs] then
                    continue;
                end if;

                J:=Matrix(k,#eqs,5,
                    &cat[[Evaluate(derivs[i][j],pt):j in [1..5]]
                          :i in [1..#eqs]]);
                r:=Rank(J);
                Jfiber:=Submatrix(J,1,2,Nrows(J),4);
                rf:=Rank(Jfiber);
                balanced+:=1;
                Append(~full_ranks,r);
                Append(~fiber_ranks,rf);

                // Pull back DB=W^2 or DC=W^2.  Since z is a unit here,
                // q is a unit times W^2 and the W=0 derivative vanishes.
                // Exact smoothness on this weighted layer is rank 3 after
                // deleting the q column.
                sr:=-1;
                if label in {"DB","DC"} then
                    sr:=rf;
                    square_balanced+:=1;
                    if sr eq 3 then square_balanced_smooth+:=1; end if;
                end if;

                de,dd:=Direction(chart,z0);
                key:=<Integers()!t0,Integers()!de,Integers()!dd,
                      Integers()!x0,Integers()!y0>;
                AddCone(~cone_keys,~cone_ranks,~cone_square_ranks,
                        key,rf,sr);
            end for; end for; end for; end for;

            printf "CHART BASE %o DEEP %o C %o RAW_BALANCED %o ",
                   label,deep,chart,balanced;
            printf "FULL_RANKS <0:%o,1:%o,2:%o,3:%o> ",
                   CountRank(full_ranks,0),CountRank(full_ranks,1),
                   CountRank(full_ranks,2),CountRank(full_ranks,3);
            printf "FIBER_RANKS <0:%o,1:%o,2:%o,3:%o> ",
                   CountRank(fiber_ranks,0),CountRank(fiber_ranks,1),
                   CountRank(fiber_ranks,2),CountRank(fiber_ranks,3);
            printf "STRICT_SMOOTH_CERTIFIED %o UNRESOLVED %o",
                   CountRank(fiber_ranks,3),
                   balanced-CountRank(fiber_ranks,3);
            if label in {"DB","DC"} then
                printf " DUAL_SQUARE_BALANCED %o DUAL_SQUARE_SMOOTH %o",
                       square_balanced,square_balanced_smooth;
            end if;
            printf "\n";
        end for;

        smooth_keys:=[cone_keys[i]:i in [1..#cone_keys]|
                      cone_ranks[i] eq 3];
        square_smooth_keys:=[];
        if label in {"DB","DC"} then
            square_smooth_keys:=[cone_keys[i]:i in [1..#cone_keys]|
                                cone_square_ranks[i] eq 3];
        end if;

        grand_cones+:=#cone_keys;
        grand_smooth+:=#smooth_keys;
        grand_square_smooth+:=#square_smooth_keys;

        printf "PAIR BASE %o DEEP %o RAW_BALANCED_CONES %o ",
               label,deep,#cone_keys;
        printf "FIBER_RANKS <0:%o,1:%o,2:%o,3:%o> STRICT_SMOOTH_CERTIFIED %o UNRESOLVED %o",
               CountRank(cone_ranks,0),CountRank(cone_ranks,1),
               CountRank(cone_ranks,2),CountRank(cone_ranks,3),
               #smooth_keys,#cone_keys-#smooth_keys;
        if label in {"DB","DC"} then
            printf " DUAL_SQUARE_SMOOTH %o",#square_smooth_keys;
        end if;
        printf "\n";
        print " SMOOTH_CONE_KEYS",smooth_keys;
        if label in {"DB","DC"} then
            print " DUAL_SQUARE_SMOOTH_CONE_KEYS",square_smooth_keys;
        end if;
    end for;
end for;

print "TOTAL_DEDUP_RAW_BALANCED_CONES",grand_cones;
print "TOTAL_STRICT_SMOOTH_CERTIFIED_BALANCED_CONES",grand_smooth;
print "TOTAL_DB_DC_DUAL_SQUARE_SMOOTH_BALANCED_CONES",grand_square_smooth;
print "CONTACT6_M612_AFFINE_FIRST_BLOWUPS_MOD5_STRICT_DONE";
quit;
