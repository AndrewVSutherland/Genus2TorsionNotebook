# p-adic boundary lift counts for H(t,r)=G(t^2,r)=0.
# These counts test whether the small-prime good-chart obstruction kills boundary lifts. It does not.

from sage.all import *

T = PolynomialRing(QQ, 'u'); u = T.gen()
K = FractionField(T)
R = PolynomialRing(K, 'x'); x = R.gen()
S = PolynomialRing(QQ, ['q','r']); q,r = S.gens()
A = PolynomialRing(ZZ, ['t','y']); t,y = A.gens()

def build_H_terms():
    beta=(u**2+1)**2/(4*u**2); sbeta=(u**2+1)/(2*u); s=(u**2-1)/(2*u)
    alpha=K(beta)-K(s**5)/K(beta*sbeta); lam=(alpha-1)**4/((alpha-K(beta))**2*alpha)
    expr=(x**3*(x-alpha)**2-(x-1)*((x-1)**4-lam*(x-K(beta))**2*x))/(2*(x-alpha)*(x-K(beta)))
    A0,rem=R(expr.numerator()).quo_rem(R(expr.denominator()))
    f=R(A0**2-lam*x**4*(x-1))
    lcmden=T(1)
    for c in f.coefficients(sparse=False): lcmden=lcm(lcmden,T(c.denominator()))
    G=S(0)
    for i,c in enumerate(f.coefficients(sparse=False)):
        poly=T(lcmden*c)
        for j,cj in enumerate(poly.list()):
            if cj: G += QQ(cj)*q**(j//2)*r**i
    lcmd=LCM([QQ(c).denominator() for c in G.coefficients()]); G=S(lcmd)*G
    cont=GCD([ZZ(c) for c in G.coefficients()]); G=G/S(cont)
    H=A(0)
    for (iq,ir), coeff in G.dict().items(): H += ZZ(coeff) * t**(2*iq) * y**ir
    return [(int(c), mon[0], mon[1]) for mon,c in H.dict().items()]

terms = build_H_terms()

def eval_terms(a,b,mod):
    s=0
    for c,i,j in terms:
        s=(s + c*pow(a,i,mod)*pow(b,j,mod)) % mod
    return s

def lift_solutions(p,N):
    sols=[]
    for a in range(p):
        for b in range(p):
            if eval_terms(a,b,p)==0:
                sols.append((a,b))
    print('p',p,'k',1,'count',len(sols),'classes',sols)
    mod=p
    for k in range(1,N):
        new=[]; mod2=mod*p
        for a,b in sols:
            for da in range(p):
                aa=a+mod*da
                for db in range(p):
                    bb=b+mod*db
                    if eval_terms(aa,bb,mod2)==0:
                        new.append((aa,bb))
        sols=new; mod=mod2
        classes={}
        for a,b in sols:
            classes[(a%p,b%p)]=classes.get((a%p,b%p),0)+1
        print('p',p,'k',k+1,'count',len(sols),'classes',classes)

lift_solutions(3,6)
lift_solutions(5,4)
