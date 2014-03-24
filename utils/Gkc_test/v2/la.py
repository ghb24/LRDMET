import numpy as N
import scipy.linalg
import math
import cmath
import utils

def peigh(a, b, thresh=1.e-12):
    """
    See related fn, scipy.linalg.eigsh.

    Solve real symmetric or complex hermitian eigenvalue problem
    when b is SINGULAR, by solving in the non-null space
    of b.
    """
    return peig(a, b, thresh, scipy.linalg.eigh)

def peig(a, b, thresh=1.e-12, scipy_eig_func=scipy.linalg.eig):
    """
    See related fn, scipy.linalg.eig

    Solve real or complex eigenvalue problem
    when b is SINGULAR, by solving in the non-null space
    of b.
    """
    beigs, bvecs=scipy_eig_func(b)
    bspace=[index for index, eig in enumerate(beigs) if abs(eig) > thresh]
    pbeigs = beigs[bspace]

    brows, bcols=b.shape

    if any([eig < 0. for eig in pbeigs]):
        pbvecs= N.empty((brows, len(bspace)), N.complex128)
        for pindex, index in enumerate(bspace):
            pbvecs[:,pindex]=bvecs[:,index]/cmath.sqrt(beigs[index])
    else:
        pbvecs= N.empty((brows, len(bspace)))
        for pindex, index in enumerate(bspace):
            pbvecs[:,pindex]=bvecs[:,index]/math.sqrt(beigs[index])

    pa=N.dot(N.dot(pbvecs.T, a), pbvecs)

    paeigs, pavecs=scipy_eig_func(pa)

    return paeigs, N.dot(pbvecs, pavecs)

def solve_perturb(omega,s,h,e0,c0,v,project=True,sign_ham=1.):

    if project:
        p0=N.outer(c0,N.conj(c0))
        p0=s-N.dot(p0,s)

    se,sv=scipy.linalg.eigh(s)
    spani=[]
    for i, sei in enumerate(se):
        if sei>1.e-10:
            spani.append(i)

    span_sv=sv[:,spani]
    if len(spani)<s.shape[0]:
        print "truncated", s.shape[0]-len(spani), "vecs"

    lhs=omega*s-sign_ham*(h-e0*s)

    # work in span of s
    lhs=N.dot(N.conj(span_sv.T), N.dot(lhs,span_sv))

    print "dimension", s.shape, "smallest eig", scipy.linalg.eigvalsh(s)[0], list(sorted(N.abs(scipy.linalg.eigvalsh(lhs))))[0]
    rhs=N.dot(v,c0)

    if project:
        #rhs-=N.dot(N.conj(c0),N.dot(s,rhs))*c0
        rhs=N.dot(p0,rhs)

    # work in span of s
    rhs=N.dot(N.conj(span_sv.T),rhs)

    sol=scipy.linalg.lstsq(lhs,rhs)
    psi1=sol[0]

    #sol=N.dot(scipy.linalg.inv(lhs),rhs)
    #psi1=sol

    # return to full space
    psi1=N.dot(span_sv,psi1)

    if project:
        #psi1=N.dot(p0,psi1)
        psi1-=N.dot(N.conj(c0),N.dot(s,psi1))*c0

        overlap=N.dot(N.conj(psi1),N.dot(s,c0))
        if abs(overlap) > 1.e-12: 
            print "warning: overlap lost", overlap, scipy.linalg.norm(psi1)
    return psi1



# def solve_perturb(omega,s,h,e0,c0,v,project=True,sign_ham=1.):

#     if project:
#         p0=N.outer(c0,N.conj(c0))
#         p0=s-N.dot(p0,s)

#     lhs=omega*s-sign_ham*(h-e0*s)

#     rhs=N.dot(v,c0)

#     if project:
#         rhs=N.dot(p0,rhs)

#     sol=scipy.linalg.lstsq(lhs,rhs)
#     psi1=sol[0]
#     if project:
#         psi1=N.dot(p0,psi1)

#     if abs(N.dot(N.conj(c0),psi1))>1.e-10: print "lost orthogonality", N.dot(N.conj(c0),psi1)
#     return psi1

def solve_perturb_complex(omega,s,h,e0,c0,v,delta=0.,sign_ham=1.):
    lhs_diag=omega*s-sign_ham*(h-e0*s)
    lhs_off=sign_ham*delta*s
    rhs_block=N.dot(v,c0)
    
    dim=s.shape[0]
    lhs_full=utils.zeros([2*dim,2*dim])
    lhs_full[0:dim,0:dim]=lhs_diag
    lhs_full[0:dim,dim:]=-lhs_off
    lhs_full[dim:,dim:]=lhs_diag
    lhs_full[dim:,0:dim]=lhs_off

    rhs_full=utils.zeros([2*dim])
    rhs_full[0:dim]=rhs_block

    #sol=scipy.linalg.lstsq(lhs_full,rhs_full)
    #psi1=sol[0]
    sol=N.dot(scipy.linalg.inv(lhs_full),rhs_full)
    psi1=sol
    return psi1[0:dim],psi1[dim:]
