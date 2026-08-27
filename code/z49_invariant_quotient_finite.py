#!/usr/bin/env python3
"""Finite geometry of the sign quotient of the order-49 incidence curve.

The four incidence equations are invariant under (u,v)->(-u,-v).
On u!=0 use w=u^2 and s=v/u.  On v!=0 use t=v^2 and q=-u/v.
Both quotient charts are derived without expansion.  Finite quotient points,
including nonsquare twists invisible in the original incidence, are enumerated
by a triangular monic-square test.  No Groebner basis or height search is used.
"""

import argparse
import resource
import sys
from itertools import product
from math import comb

import sympy as sp

sys.path.insert(0,"code")
import z49_structural_contact_iterate as contact


def quotient_polynomials(equations,chart):
    a,b,u,v,r=sp.symbols("a b u v r")
    if chart=="slope":
        w,s=sp.symbols("w s")
        variables=(a,b,w,s,r)
    else:
        t,q=sp.symbols("t q")
        variables=(a,b,t,q,r)
    out=[]
    for source in equations:
        expr=0
        for monomial,coefficient in source.terms():
            ia,ib,iu,iv,ir=monomial
            assert (iu+iv)%2==0
            if chart=="slope":
                expr += coefficient*a**ia*b**ib*w**((iu+iv)//2)*s**iv*r**ir
            else:
                expr += coefficient*((-1)**iu)*a**ia*b**ib*t**((iu+iv)//2)*q**iu*r**ir
        out.append(sp.Poly(expr,*variables,domain=sp.ZZ))
    return variables,out


def trim(poly,p):
    ans=[c%p for c in poly]
    while len(ans)>1 and ans[-1]==0:
        ans.pop()
    return ans


def gcd_poly(a,b,p):
    a,b=trim(a,p),trim(b,p)
    while b!=[0]:
        rr=a[:]
        while rr!=[0] and len(rr)>=len(b):
            shift=len(rr)-len(b)
            scale=rr[-1]*pow(b[-1],-1,p)%p
            for i,c in enumerate(b):
                rr[i+shift]=(rr[i+shift]-scale*c)%p
            rr=trim(rr,p)
        a,b=b,rr
    return a


def f_coefficients(a,b,p):
    inv4=pow(4,-1,p)
    return [(2*a-35*inv4)%p,
            (2*b-7*a+35)%p,
            (a*a-7*b-35)%p,
            (2*a*b+21)%p,
            (b*b-7)%p,1]


def smooth_f(f,p):
    derivative=[i*f[i]%p for i in range(1,len(f))]
    return len(gcd_poly(f,derivative,p))==1


def multiply(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        for j,bj in enumerate(b):
            out[i+j]=(out[i+j]+ai*bj)%p
    return out


def add_polys(a,b,p):
    n=max(len(a),len(b)); out=[0]*n
    for i in range(n):
        out[i]=((a[i] if i<len(a) else 0)+
                (b[i] if i<len(b) else 0))%p
    return out


def rhs_coefficients(r,p):
    base=[comb(7,k)*pow(-r,7-k,p)%p for k in range(8)]
    out=[0]*9
    for k,c in enumerate(base):
        out[k]=(out[k]-c)%p
        out[k+1]=(out[k+1]+c)%p
    return out


def is_monic_quartic_square(g,p):
    """Test g=A^2, A monic quartic, by triangular square rooting."""
    if len(g)!=9 or g[8]%p!=1:
        return False
    inv2=pow(2,-1,p)
    z3=g[7]*inv2%p
    z2=(g[6]-z3*z3)*inv2%p
    z1=(g[5]-2*z3*z2)*inv2%p
    z0=(g[4]-2*z3*z1-z2*z2)*inv2%p
    aa=[z0,z1,z2,z3,1]
    return multiply(aa,aa,p)==[c%p for c in g]


def enumerate_quotient(p):
    """Return every open affine quotient point, including nonsquare twists."""
    inv2=pow(2,-1,p)
    curves=[]
    for a,b in product(range(p),repeat=2):
        if (a+b-5*inv2)%p==0:
            continue
        f=f_coefficients(a,b,p)
        if smooth_f(f,p):
            curves.append((a,b,f))
    rhs={r:rhs_coefficients(r,p) for r in range(p) if r!=1}
    points=[]

    # v != 0: B^2=t*(x-q)^2.  The invariant t need not be a square.
    for a,b,f in curves:
        for q in range(p):
            if q==1:
                continue
            base=multiply([(q*q)%p,(-2*q)%p,1],f,p)
            for t in range(1,p):
                tf=[t*c%p for c in base]
                for r,rr in rhs.items():
                    if r==q:
                        continue
                    if is_monic_quartic_square(add_polys(tf,rr,p),p):
                        points.append(("root",a,b,t,q,r))

    # v=0,u!=0: B^2=w; this is the s=0 chart.
    for a,b,f in curves:
        for w in range(1,p):
            wf=[w*c%p for c in f]
            for r,rr in rhs.items():
                if is_monic_quartic_square(add_polys(wf,rr,p),p):
                    points.append(("constant",a,b,w,0,r))

    assert len(points)==len(set(points))
    return points


def quotient_projections(points,p):
    projections={"qt":set(),"sb":set(),"sr":set()}
    liftable=0
    for chart,a,b,z,c,r in points:
        if chart=="root":
            t,q=z,c
            projections["qt"].add((q,t))
            if pow(t,(p-1)//2,p)==1:
                liftable += 2
            if q!=0:
                s=(-pow(q,-1,p))%p
                projections["sb"].add((s,b))
                projections["sr"].add((s,r))
        else:
            w=z
            if pow(w,(p-1)//2,p)==1:
                liftable += 2
            projections["sb"].add((0,b))
            projections["sr"].add((0,r))
    return projections,liftable


def rank_mod(rows,p):
    if not rows:
        return 0
    mat=[row[:] for row in rows]
    nr,nc=len(mat),len(mat[0]); rank=0
    for col in range(nc):
        pivot=next((i for i in range(rank,nr) if mat[i][col]%p),None)
        if pivot is None:
            continue
        mat[rank],mat[pivot]=mat[pivot],mat[rank]
        inv=pow(mat[rank][col],-1,p)
        mat[rank]=[x*inv%p for x in mat[rank]]
        for i in range(nr):
            if i==rank or mat[i][col]%p==0:
                continue
            scale=mat[i][col]%p
            mat[i]=[(mat[i][j]-scale*mat[rank][j])%p
                    for j in range(nc)]
        rank+=1
        if rank==nr:
            break
    return rank


def monomials_total(degree):
    return [(i,j) for i in range(degree+1)
            for j in range(degree+1-i)]


def monomials_box(dx,dy):
    return [(i,j) for i in range(dx+1) for j in range(dy+1)]


def evaluation_rank(points,monomials,p):
    rows=[[pow(x,i,p)*pow(y,j,p)%p for i,j in monomials]
          for x,y in sorted(points)]
    return rank_mod(rows,p)


EXPECTED_LIFT_COUNTS={5:4,11:10,13:12,17:18}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--memory-mb",type=int,default=180)
    parser.add_argument("--primes",default="5,11,13,17")
    parser.add_argument("--max-total-degree",type=int,default=4)
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,
                       (args.memory_mb*1024**2,args.memory_mb*1024**2))
    _x,_variables,_h,_f,_solved,equations=contact.derive_system()
    print("Z49_INVARIANT_QUOTIENT_GEOMETRY")
    for chart in ("slope","root"):
        variables,polys=quotient_polynomials(equations,chart)
        print("chart",chart,"variables",variables)
        for i,poly in enumerate(polys):
            print(" equation",i,"degrees",poly.degree_list(),
                  "total_degree",poly.total_degree(),"terms",len(poly.terms()))

    prime_data={}
    for p in [int(x) for x in args.primes.split(",") if x.strip()]:
        points=enumerate_quotient(p)
        projections,liftable=quotient_projections(points,p)
        if p in EXPECTED_LIFT_COUNTS:
            assert liftable==EXPECTED_LIFT_COUNTS[p]
        root_count=sum(pt[0]=="root" for pt in points)
        constant_count=len(points)-root_count
        print("QUOTIENT",p,"points",len(points),"root",root_count,
              "constant",constant_count,"lifted_open_incidence",liftable)
        print(" samples",points[:8])
        for name in ("qt","sb","sr"):
            print(" projection",name,"distinct",len(projections[name]),
                  "points",sorted(projections[name]))
        prime_data[p]=projections

    for name in ("qt","sb","sr"):
        print("RANK_DIAGNOSTIC",name)
        for degree in range(1,args.max_total_degree+1):
            mons=monomials_total(degree)
            rows=[]
            for p,projections in prime_data.items():
                rank=evaluation_rank(projections[name],mons,p)
                rows.append((p,len(projections[name]),rank,len(mons)-rank))
            print(" total_degree",degree,"monomials",len(mons),rows)
        for dx,dy in ((1,1),(1,2),(2,1),(2,2),(1,3),(3,1)):
            mons=monomials_box(dx,dy)
            rows=[]
            for p,projections in prime_data.items():
                rank=evaluation_rank(projections[name],mons,p)
                rows.append((p,len(projections[name]),rank,len(mons)-rank))
            print(" bidegree",(dx,dy),"monomials",len(mons),rows)
    print("Z49_INVARIANT_QUOTIENT_GEOMETRY_DONE")


if __name__=="__main__":
    main()
