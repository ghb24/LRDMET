#makes configs
import numpy as N
import numpy.linalg
import models
import fermions
import utils

def _spinify(mat):
    sp_mat=utils.zeros([mat.shape[0]*2,mat.shape[1]*2])
    for i in xrange(mat.shape[0]):
        for j in xrange(mat.shape[1]):
            sp_mat[2*i,2*j]=mat[i,j]            
            sp_mat[2*i+1,2*j+1]=mat[i,j]
    return sp_mat

def _spinify_vec(vec):
    sp_vec=utils.zeros([vec.shape[0]*2])
    for i in xrange(vec.shape[0]):
        sp_vec[2*i]=vec[i]            
        #sp_vec[2*i+1]=vec[i]
    return sp_vec
    
def _add_core(configs, core):
    new_configs=[config+core for config in configs]
    return new_configs

def _unspinify(sp_mat):
    mat=utils.zeros([sp_mat.shape[0]/2,sp_mat.shape[1]/2])
    for i in xrange(mat.shape[0]):
        for j in xrange(mat.shape[1]):
            mat[i,j]=sp_mat[2*i,2*j]
    return mat
    
def _linear_orthonormal_basis(a, orient='c',thresh=1.e-12):
    """
    create linearly independent basis from columns
    (optionally rows) of a
    """
    if orient=='r':
        inp_a=a.T
    else:
        inp_a=a
    s=N.dot(N.conj(inp_a.T),inp_a)
    eigs,vecs=numpy.linalg.eigh(s)
    keep_indices=[]
    for i, eig in enumerate(eigs):
        if abs(eig)>thresh:
            keep_indices.append(i)
    keep_vecs=vecs[:,keep_indices]
    ret_a=N.dot(inp_a, keep_vecs)
    if orient=='r':
        ret_a=ret_a.T

    return ret_a

def make_t1amp_ao(h, omega, perturb, nocc):
    eigs,cmo=numpy.linalg.eigh(h)
    nsites=h.shape[0]

    perturb_mo=N.dot(cmo.T,N.dot(perturb,cmo))

    t1amp=utils.zeros([nsites,nsites])
    for i in xrange(nocc):
        for a in xrange(nocc,nsites):
            t1amp[a,i]=perturb_mo[i,a]/(omega-(eigs[a]-eigs[i]))

    t1amp_ao=N.dot(cmo,N.dot(t1amp,cmo.T))
    return t1amp_ao

def make_greens_ao(h, omega, perturb, nocc):
    eigs,cmo=numpy.linalg.eigh(h)
    nsites=h.shape[0]

    perturb_mo=N.dot(cmo.T,perturb)
    #damp=utils.zeros([nocc])
    #camp=utils.zeros([nsites-nocc])
    
    damp=utils.zeros([nsites])
    camp=utils.zeros([nsites])

    for i in xrange(nocc):
        damp[i]=perturb_mo[i]/(omega-eigs[i])
        #print i, perturb_mo[i]/(omega-eigs[i]), damp[i]
    for a in xrange(nocc,nsites):
        camp[a]=perturb_mo[a]/(omega-eigs[a])

    #damp_ao=N.dot(cmo[:,:nocc],damp)
    damp_ao=N.dot(cmo,damp)
    #camp_ao=N.dot(cmo[:,nocc:],camp)
    camp_ao=N.dot(cmo,camp)

    return damp_ao,camp_ao

def mb_sp_ops(hlat,perturb,omega,nimp,nocc,pemb,pcore,pvirt):
    """
    single particle external ops 
    (for single particle Green's function)
    """
    nimp_sp=nimp*2
    nsites_sp=hlat.shape[0]*2
    nocc_sp=nocc*2
    
    # transform h into imp+bath/ext(all) basis
    pall=N.hstack([pemb,pcore,pvirt])
    hall=N.dot(pall.T,N.dot(hlat,pall))
    hall_sp=_spinify(hall)
    
    # green's function amps in (all) basis
    damp_ao,camp_ao=make_greens_ao(hlat,omega,perturb,nocc)
    damp_all=N.dot(pall.T,damp_ao)
    camp_all=N.dot(pall.T,camp_ao)
    #print N.dot(damp_ao,pall)
    damp_all_sp=_spinify_vec(damp_all)
    camp_all_sp=_spinify_vec(camp_all)

    # mb_ext operators
    sites_imp=range(0,nimp_sp*2)
    sites_ext=range(nimp_sp*2,nsites_sp)
    
    coeff_c=camp_all_sp[sites_ext]
    cops=[models.ContractedC(sites_ext,coeff_c)]

    coeff_d=damp_all_sp[sites_ext]
    dops=[models.ContractedD(sites_ext,coeff_d)]

    ops={ "1":models.Unit(), 
          "c":cops, 
          "d":dops }

    return ops

def mb_ph_ops(hlat,perturb,omega,nimp,nocc,pemb,pcore,pvirt):
    """
    external ops for part-hole Green's function
    """
    nimp_sp=nimp*2
    nsites_sp=hlat.shape[0]*2
    nocc_sp=nocc*2

    # transform h into imp+bath/ext(all) basis
    pall=N.hstack([pemb,pcore,pvirt])
    hall=N.dot(pall.T,N.dot(hlat,pall))
    hall_sp=_spinify(hall)

    # t1amp in (all) basis
    t1amp_ao=make_t1amp_ao(hlat,omega,perturb,nocc)
    t1amp_all=N.dot(pall.T,N.dot(t1amp_ao,pall))
    t1amp_all_sp=_spinify(t1amp_all)

    # mb_ext_operators
    sites_imp=range(0,nimp_sp*2)
    sites_ext=range(nimp_sp*2,nsites_sp)
    
    cops=[]
    dops=[]

    project=False

    if project:
        li_t1amp_ext_imp=_linear_orthonormal_basis(t1amp_all_sp[sites_ext,:][:,sites_imp],orient='c')
        print "retained %i out of %i vectors" % (li_t1amp_ext_imp.shape[1],len(sites_imp))

        li_t1amp_imp_ext=_linear_orthonormal_basis(t1amp_all_sp[sites_imp,:][:,sites_ext],orient='r')
        print "retained %i out of %i vectors" % (li_t1amp_imp_ext.shape[0],len(sites_imp))
                                                 
        for i in xrange(li_t1amp_ext_imp.shape[1]):
            cops.append(models.ContractedC(sites_ext,li_t1amp_ext_imp[:,i]))
            
        for i in xrange(li_t1amp_imp_ext.shape[0]):
            dops.append(models.ContractedD(sites_ext,li_t1amp_imp_ext[i,:]))
    else:
        for imp in sites_imp:
            coeff_c=t1amp_all_sp[sites_ext,imp]
            cops.append(models.ContractedC(sites_ext,coeff_c))
            coeff_d=t1amp_all_sp[imp,sites_ext]
            dops.append(models.ContractedD(sites_ext,coeff_d))

    cd_coeffs={}
    for i in sites_ext:
        for j in sites_ext:
            if abs(t1amp_all_sp[i,j]) > 1.e-8:
                cd_coeffs[i,j]=t1amp_all_sp[i,j]
    cdop=models.ContractedCD(cd_coeffs)

    ops={ "1":models.Unit(), 
          "c":cops, 
          "d":dops, 
          "cd":cdop }

    return ops

def mb_configs(nsites,nimp,ne_imp_sp,ncore_sp,sz):
    nsites_sp=nsites*2
    ndim=2
    nimp_sp=nimp*2

    # fermion spaces
    n_configs=fermions.conserving_configs(nimp_sp*2,ndim,ne_imp_sp,sz)
    np_configs=fermions.conserving_configs(nimp_sp*2,ndim,ne_imp_sp+1,sz+1)
    np_configs+=fermions.conserving_configs(nimp_sp*2,ndim,ne_imp_sp+1,sz-1)
    nm_configs=fermions.conserving_configs(nimp_sp*2,ndim,ne_imp_sp-1,sz+1)
    nm_configs+=fermions.conserving_configs(nimp_sp*2,ndim,ne_imp_sp-1,sz-1)

    all_configs=n_configs+np_configs+nm_configs

    configs_dict={ "n":n_configs, 
                   "n+1":np_configs, 
                   "n-1":nm_configs,
                   "all":all_configs }

    # add core determinant
    nvirt_sp=nsites_sp-2*nimp_sp-ncore_sp
    core_config=[1]*ncore_sp+[0]*nvirt_sp
    for key in configs_dict:
        configs_dict[key]=_add_core(configs_dict[key], core_config)

    return configs_dict

def mb_ph_ops_configs(ops_dict,configs_dict):
    ops_configs=[(ops_dict["1"],configs_dict["n"])]
    ops_configs+=[(ops_dict["cd"],configs_dict["n"])]

    for op in ops_dict["c"]:
        ops_configs+=[(op,configs_dict["n-1"])]

    for op in ops_dict["d"]:
        ops_configs+=[(op,configs_dict["n+1"])]
    return ops_configs

#creates neutral ops configs (N particles - operator is unity)
def mb_sp_ops_configs(ops_dict,configs_dict):
    neutral_op_configs=[(ops_dict["1"],configs_dict["n"])]

    plus_op_configs=[(ops_dict["1"],configs_dict["n+1"])]
    for op in ops_dict["c"]:
        plus_op_configs+=[(op,configs_dict["n"])]
    
    minus_op_configs=[(ops_dict["1"],configs_dict["n-1"])]
    for op in ops_dict["d"]:
        minus_op_configs+=[(op,configs_dict["n"])]

    return neutral_op_configs,plus_op_configs,minus_op_configs
