import numpy as N
import embed
import models
import utils

def _spinify(mat):
    sp_mat=utils.zeros([mat.shape[0]*2,mat.shape[1]*2])
    for i in xrange(mat.shape[0]):
        for j in xrange(mat.shape[1]):
            sp_mat[2*i,2*j]=mat[i,j]            
            sp_mat[2*i+1,2*j+1]=mat[i,j]
    return sp_mat

# half-filled Hubbard 1-site DMET model
def hubbard_imp_ext_ni_h(t,u,nocc,nimp,nsites):
    ncore=2*nocc-2*nimp
    cocc=range(4*nimp,4*nimp+ncore)

    hlat=utils.zeros([nsites,nsites])
    for i in xrange(nsites):
        hlat[i,i]=+u/2
        for j in xrange(nsites):
            if abs(i-j)==1:
                hlat[i,j]=t
                
    hlat_sp=_spinify(hlat)    

    pemb,pcore,pvirt=embed.embed(hlat,nimp,nocc)
    
    pall=N.hstack([pemb,pcore,pvirt])
    hall=N.dot(pall.T,N.dot(hlat,pall))
    hall_sp=_spinify(hall)

    hall_sp_diag=N.diag(hall_sp)
    e0=sum(hall_sp_diag[cocc])
    return hlat, hall, hlat_sp, hall_sp, e0


def hubbard_imp_ext_h(h0,u,nocc,nimp,nsites):
    nsites_sp=nsites*2
    sites_imp=range(0,4*nimp)
    sites_ext=range(4*nimp,2*nsites)

    t_dict={}
    for i in xrange(nsites_sp):
        for j in xrange(nsites_sp):
            if abs(h0[i,j]) > 1.e-12:
                t_dict[i,j]=h0[i,j]
    
    t_dict[0,0]=0
    t_dict[1,1]=0

    u_dict={}
    u_dict[0,1]=u

    # h in imp+ext space
    hop=models.GeneralHubbard(t_dict,u_dict)

    t_imp={}
    for i in sites_imp:
        for j in sites_imp:
            if abs(h0[i,j]) > 1.e-12:
                t_imp[i,j]=h0[i,j]

    t_imp[0,0]=0
    t_imp[1,1]=0

    t_ext={}
    for i in sites_ext:
        for j in sites_ext:
            if abs(h0[i,j]) > 1.e-12:
                t_ext[i,j]=h0[i,j]

    hext=models.ContractedCD(t_ext)

    himp=models.GeneralHubbard(t_imp,u_dict)

    hcs_imp=[]
    hds_imp=[]
    hcs_ext=[]
    hds_ext=[]
    
    for i in sites_imp:
        imp_coeffs=utils.zeros(len(sites_imp))
        imp_coeffs[i]=1.
        hcs_imp.append(models.ContractedC(sites_imp,imp_coeffs))
        hds_imp.append(models.ContractedD(sites_imp,imp_coeffs))
        hcs_ext.append(models.ContractedC(sites_ext,h0[sites_ext,i]))
        hds_ext.append(models.ContractedD(sites_ext,h0[i,sites_ext]))

    return hop, himp, hcs_imp, hds_imp, hcs_ext, hds_ext, hext



# half-filled single-site anderson
def si_anderson_imp_ext_ni_h(t,nocc,nimp,nsites):
    ncore=2*nocc-2*nimp
    cocc=range(4*nimp,4*nimp+ncore)

    hlat=utils.zeros([nsites,nsites])
    for i in xrange(nsites):
        for j in xrange(nsites):
            if abs(i-j)==1:
                hlat[i,j]=t

    hlat_sp=_spinify(hlat)    

    pemb,pcore,pvirt=embed.embed(hlat,nimp,nocc)
    
    pall=N.hstack([pemb,pcore,pvirt])
    hall=N.dot(pall.T,N.dot(hlat,pall))
    hall_sp=_spinify(hall)

    hall_sp_diag=N.diag(hall_sp)
    e0=sum(hall_sp_diag[cocc])

    return hlat, hall, hlat_sp, hall_sp, e0

def si_anderson_imp_ext_h(h0,u,nocc,nimp,nsites):
    nsites_sp=nsites*2
    sites_imp=range(0,4*nimp)
    sites_ext=range(4*nimp,2*nsites)

    t_dict={}
    for i in xrange(nsites_sp):
        for j in xrange(nsites_sp):
            if abs(h0[i,j]) > 1.e-12:
                t_dict[i,j]=h0[i,j]
    
    t_dict[0,0]=-u/2
    t_dict[1,1]=-u/2

    u_dict={}
    u_dict[0,1]=u

    # h in imp+ext space
    hop=models.GeneralHubbard(t_dict,u_dict)

    t_imp={}
    for i in sites_imp:
        for j in sites_imp:
            if abs(h0[i,j]) > 1.e-12:
                t_imp[i,j]=h0[i,j]

    t_imp[0,0]=-u/2
    t_imp[1,1]=-u/2

    t_ext={}
    for i in sites_ext:
        for j in sites_ext:
            if abs(h0[i,j]) > 1.e-12:
                t_ext[i,j]=h0[i,j]

    hext=models.ContractedCD(t_ext)

    himp=models.GeneralHubbard(t_imp,u_dict)

    hcs_imp=[]
    hds_imp=[]
    hcs_ext=[]
    hds_ext=[]
    
    for i in sites_imp:
        imp_coeffs=utils.zeros(len(sites_imp))
        imp_coeffs[i]=1.
        hcs_imp.append(models.ContractedC(sites_imp,imp_coeffs))
        hds_imp.append(models.ContractedD(sites_imp,imp_coeffs))
        hcs_ext.append(models.ContractedC(sites_ext,h0[sites_ext,i]))
        hds_ext.append(models.ContractedD(sites_ext,h0[i,sites_ext]))

    return hop, himp, hcs_imp, hds_imp, hcs_ext, hds_ext, hext

