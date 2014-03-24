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

    lhs=omega*s-sign_ham*(h-e0*s)

    rhs=N.dot(v,c0)

    if project:
        rhs=N.dot(p0,rhs)

    sol=scipy.linalg.lstsq(lhs,rhs)
    psi1=sol[0]
    if project:
        psi1=N.dot(p0,psi1)

    return psi1

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
