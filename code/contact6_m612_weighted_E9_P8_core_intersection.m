//////////////////////////////////////////////////////////////////////
// Bounded exact fiber-product search: rational endpoint R3 P8 family
// against the independent cubic-contact [3,3] core on b=0, a=1/e.
//
// The P8 normalization is parametrized by tau.  We deduplicate its e-values,
// solve the fully saturated fixed-curve core exactly, require M=L^2, and
// verify the two independent 3-directions in the Jacobian.  Any survivor is
// then audited on both Richelot twists with the exact dual-class test.
//
// Usage:
//   magma -b height:=30 code/contact6_m612_weighted_E9_P8_core_intersection.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned height then height:=30;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned progress_interval then progress_interval:=50;
elif Type(progress_interval) eq MonStgElt then
    progress_interval:=StringToInteger(progress_interval);
end if;

Q:=Rationals(); Z:=Integers();
load "code/contact6_m612_tb_core_tools.m";
run_known_seed:=false; NoMain:=true;
load "code/contact6_m612_dual_class_exact.m";

// Exact rational parametrization of the full P8 space component.
K<tau>:=FunctionField(Q);
t:=4*(tau^2+tau-6)/(tau^2+6);
y:=2*(tau^2-24*tau-6)/(tau^2+6);
assert y^2 eq 100-6*t^2;
den:=t^8-50*t^6+Q!4375/3*t^4-Q!62500/3*t^2+Q!1562500/9;
eP:=-(Q!25/3)*t^2/(t^4-25*t^2+Q!1250/3);
muP:=(-5*t^7+Q!1750/9*t^5-Q!6250/3*t^3)/den;
Bnu:=(-Q!500/3*t^5+Q!175000/27*t^3-Q!625000/9*t)/den;
hdisc:=(t^7-Q!500/9*t^5+Q!28750/27*t^3-Q!62500/9*t)/den;
nuP:=(-Bnu+hdisc*y)/2;

function ParametersOfHeight(B)
    vals:=[]; seen:={};
    for d in [1..B] do for n in [-B..B] do
        if GCD(n,d) ne 1 then continue; end if;
        q:=Q!n/d; key:=Sprint(q);
        if key in seen then continue; end if;
        Include(~seen,key); Append(~vals,q);
    end for; end for;
    return vals;
end function;

function EvaluateFunction(q,x)
    num:=Evaluate(Numerator(q),x); denq:=Evaluate(Denominator(q),x);
    if denq eq 0 then return false,Q!0; end if;
    return true,Q!num/Q!denq;
end function;

function EvaluateFunctionAtInfinity(q)
    num:=Numerator(q); denq:=Denominator(q);
    dn:=Degree(num); dd:=Degree(denq);
    if dn gt dd then return false,Q!0; end if;
    if dn lt dd then return true,Q!0; end if;
    return true,Q!LeadingCoefficient(num)/Q!LeadingCoefficient(denq);
end function;

params:=ParametersOfHeight(height);
seen_e:={}; tested:=0; duplicate_e:=0; parameter_poles:=0; e_zero:=0;
pole_data:=[]; e_zero_data:=[]; exact_empty:=0;
sat_fail:=0; variety_fail:=0; exceptional_fibers:=[];
core_q_points:=0; square_M:=0; verified_core:=0; dual_audits:=0;
quotient_only:=[]; hits:=[];
start:=Cputime();

print "CONTACT6_M612_WEIGHTED_E9_P8_CORE_INTERSECTION";
print "HEIGHT",height,"PARAMETER_COUNT",#params;
print "P8_PARAMETER","t=4(tau^2+tau-6)/(tau^2+6)";

ok_einf,einf:=EvaluateFunctionAtInfinity(eP);
ok_minf,muinf:=EvaluateFunctionAtInfinity(muP);
ok_ninf,nuinf:=EvaluateFunctionAtInfinity(nuP);
assert ok_einf and ok_minf and ok_ninf;
print "PROJECTIVE_TAU_INFINITY_P8_POINT",<einf,muinf,nuinf>;

for tau0 in params do
    oke,e0:=EvaluateFunction(eP,tau0);
    okm,mu0:=EvaluateFunction(muP,tau0);
    okn,nu0:=EvaluateFunction(nuP,tau0);
    if not (oke and okm and okn) then
        parameter_poles+:=1; Append(~pole_data,tau0); continue;
    end if;
    if e0 eq 0 then
        e_zero+:=1; Append(~e_zero_data,<tau0,e0,mu0,nu0>); continue;
    end if;
    key:=Sprint(e0);
    if key in seen_e then duplicate_e+:=1; continue; end if;
    Include(~seen_e,key); tested+:=1;
    a0:=1/e0; b0:=Q!0;
    sols,dim,satok,varok:=M612_FixedCoreSolutions(a0,b0);
    if not satok then sat_fail+:=1; end if;
    if dim eq -1 then
        exact_empty+:=1;
        continue;
    elif dim gt 0 then
        Append(~exceptional_fibers,<tau0,e0,dim,satok,varok>);
        continue;
    elif dim ne 0 then
        Append(~exceptional_fibers,<tau0,e0,dim,satok,varok>);
        continue;
    end if;
    if not varok then
        variety_fail+:=1;
        continue;
    end if;
    core_q_points+:=#sols;
    for pt in sols do
        M0:=Q!pt[1]; U0:=Q!pt[2]; v0:=Q!pt[3];
        sq,L0:=IsSquare(M0);
        if not sq or M0 eq 0 then
            if #quotient_only lt 40 then
                Append(~quotient_only,<tau0,e0,mu0,nu0,M0,U0,v0>);
            end if;
            continue;
        end if;
        square_M+:=1;
        okcore,L,ordD,ordE,Eclass:=M612_VerifyCorePoint(a0,b0,M0,U0,v0);
        if not okcore then continue; end if;
        verified_core+:=1;
        dual_records:=[];
        for eps in [1,-1] do
            okdual,rec:=ExactDualClassRecord(a0,b0,eps);
            Append(~dual_records,<eps,okdual,rec>); dual_audits+:=1;
        end for;
        Append(~hits,<tau0,e0,mu0,nu0,M0,L,U0,v0,ordD,ordE,
                      dual_records>);
        print "VERIFIED_FIBER_PRODUCT_POINT",hits[#hits];
    end for;
    if progress_interval gt 0 and tested mod progress_interval eq 0 then
        print "PROGRESS","tested_e",tested,"tau",tau0,"e",e0,
              "core_Q",core_q_points,"square_M",square_M,
              "verified",verified_core,"seconds",Cputime(start);
    end if;
end for;

// The known P8 point tau=12 has a universal degenerate core section only:
// (M,U,v)=(1,-2,1), for which U^2-4v^2=0.  The saturated solver removes it.
tau12:=Q!12;
_,e12:=EvaluateFunction(eP,tau12);
_,mu12:=EvaluateFunction(muP,tau12);
_,nu12:=EvaluateFunction(nuP,tau12);
assert <e12,mu12,nu12> eq
       <Q!-200/409,Q!-36320/167281,Q!38136/167281>;
assert einf eq e12;

// Audit the known fiber before the contact-open saturation.  Its only
// rational M != 0 point is the universal repeated-root section.
RF<Md,Ud,vd>:=PolynomialRing(Q,3,"grevlex");
Nd:=(3*Ud+6*Md)/2;
Rd:=(3*Ud^2+3*vd^2+(2/e12-15)*Md-Nd^2)/2;
raw_eqs:=[2*vd^3+2*Nd*Rd-Ud^3-6*Ud*vd^2-22*Md,
          Rd^2+2*Nd*vd^3-3*Ud^2*vd^2-3*vd^4-(1/e12^2-15)*Md,
          2*Rd*vd^3-3*Ud*vd^4-(2/e12+6)*Md];
Iraw12:=Saturation(ideal<RF|raw_eqs>,ideal<RF|Md>);
raw12:=Variety(Iraw12);
assert raw12 eq [<Q!1,Q!-2,Q!1>];

print "DONE";
print "UNIQUE_E_TESTED",tested,"DUPLICATE_E",duplicate_e,
      "PARAMETER_POLES",parameter_poles,"E_ZERO",e_zero;
print "EXACT_EMPTY_OPEN_CORE_FIBERS",exact_empty;
print "SATURATION_FAILURES",sat_fail,"VARIETY_FAILURES",variety_fail,
      "EXCEPTIONAL_FIBERS",#exceptional_fibers;
print "CORE_QUOTIENT_Q_POINTS",core_q_points,"SQUARE_M",square_M,
      "VERIFIED_INDEPENDENT_CORE",verified_core,
      "DUAL_AUDITS",dual_audits,"HITS",#hits;
print "EXCEPTIONAL_FIBER_DATA",exceptional_fibers;
print "PARAMETER_POLE_DATA",pole_data;
print "E_ZERO_BOUNDARY_DATA",e_zero_data;
print "QUOTIENT_ONLY_SAMPLES",quotient_only;
print "TAU12_P8_POINT",<e12,mu12,nu12>;
print "TAU_INFINITY_DUPLICATE_E",einf eq e12,
      "TAU_INFINITY_P8_POINT",<einf,muinf,nuinf>;
print "TAU12_DEGENERATE_CORE_SECTION",<Q!1,Q!-2,Q!1>,
      "DISCRIMINANT",Q!((-2)^2-4);
print "SECONDS",Cputime(start);
print "CONTACT6_M612_WEIGHTED_E9_P8_CORE_INTERSECTION_DONE";
quit;
