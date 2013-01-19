import numpy as N
import numpy.linalg
import utils

def pimp(nimp, nlat):
    ret=utils.zeros([nlat,nlat])
    for i in xrange(nimp):
        ret[i,i]=1.
    return ret

def embed(hlat,nimp,nocc):
    nlat=hlat.shape[0]
    elat, clat=numpy.linalg.eigh(hlat)
    cocc=clat[:,0:nocc]
    dimp=N.dot(cocc,cocc.T)[:nimp,:nimp]

    # projected overlap
    # Ct P C
    s_imp=N.dot(cocc.T,N.dot(pimp(nimp,nlat),cocc))
    sigma, u=numpy.linalg.eigh(s_imp)

    # rotated basis
    cbar=N.dot(cocc,u)

    # embedding projector pemb
    pemb=utils.zeros([nlat,nimp*2])

    for i in xrange(nimp):
        pemb[i,i]=1.
    cbath=N.dot(N.eye(nlat)-pimp(nimp,nlat),cbar)

    iptr=0
    core=[]
    for i, s in enumerate(sigma):
        # entangled orbitals
        if abs(s)>1.e-8:
            pemb[:,nimp+iptr]=cbath[:,i]/(1-s)**0.5
            iptr+=1
        else:
            core.append(i)

    # core projector pcore
    pcore=cbath[:,core]

    # virtual projector pvirt
    # generate from the null space of pemb+pcore
    pemb_core=N.hstack([pemb,pcore])
    dm_emb_core=N.dot(pemb_core, pemb_core.T)
    wts,vecs=numpy.linalg.eigh(dm_emb_core)

    pvirt=utils.zeros([nlat,nlat-nocc-nimp])
    iptr=0
    for i, w in enumerate(wts):
        if abs(w)<1.e-8:
            pvirt[:,iptr]=vecs[:,i]
            iptr+=1
            
    return pemb, pcore, pvirt
