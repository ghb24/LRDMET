import math
import time
import numpy as N
import scipy.linalg
import models
import models_embed
import embed
import qoperator
import response_embed
import opbasis
import opbasis_ni
import la
import utils

def _sp_gf(h0,omega,perturb,nocc):
    # noninteracting single-particle g.f.
    # \sum_i 1/(\omega-\epsilon_i)
    eigs,cmo=scipy.linalg.eigh(h0)
    nsites=h0.shape[0]

    perturb_mo=N.dot(cmo.T,perturb)

    d_gf=0.
    c_gf=0.
    for i in xrange(nocc):
        d_gf+=perturb_mo[i]*perturb_mo[i]/(omega-eigs[i])
    for i in xrange(nocc,nsites):
        c_gf+=perturb_mo[i]*perturb_mo[i]/(omega-eigs[i])
    return d_gf,c_gf

def _second_order_energy(h0,nocc,perturb,omega):
    # non-interacting density-density response function
    # perturb is a 1-particle [nsites,nsites] perturbation matrix
    # if perturb[0,0]=1 and all else are 0, this corresponds to
    # density-density response
    eigs,cmo=scipy.linalg.eigh(h0)
    nsites=h0.shape[0]

    perturb_mo=N.dot(cmo.T,N.dot(perturb,cmo))
    t1amp=utils.zeros([nsites,nsites])
    ret_val=0.
    for i in xrange(nocc):
        for a in xrange(nocc,nsites):
            t1amp[a,i]=perturb_mo[a,i]/(omega-(eigs[a]-eigs[i]))
            ret_val+=t1amp[a,i]*perturb_mo[i,a]
    return ret_val

def ph_greens():
    # main loop for general density-density (ph) response functions
    utils.dtype=N.complex128

    nimp=1
    nimp_sp=2*nimp
    nocc=20
    nsites=40
    nsites_sp=nsites*2
    ndim=2
    sz=0
    # sites numbered in order as imp+bath ... ext
    sites_imp=range(0,2*nimp_sp)
    sites_ext=range(2*nimp_sp,nsites_sp)
    ncore=2*nocc-2*nimp
    cocc=range(2*nimp_sp,2*nimp_sp+ncore)
    vocc=range(2*nimp_sp+ncore,nsites_sp)

    N.set_printoptions(precision=3)

    t=-1. # hopping
    delta=0.25 # broadening

#    for u in [0.0,4.0,10.0]:
    for u in [0.0,1.0,4.0,8.0]:

        # Single impurity Anderson model
        mu=0.
        hlat,hall,hlat_sp,hall_sp,e0=models_embed.si_anderson_imp_ext_ni_h(t,nocc,nimp,nsites)
        hop,himp,hcs_imp,hds_imp,hcs_ext,hds_ext,hext=models_embed.si_anderson_imp_ext_h(hall_sp,u,nocc,nimp,nsites)

#        single site Hubbard in DMET basis
#        mu=u/2
#        hlat,hall,hlat_sp,hall_sp,e0=models_embed.hubbard_imp_ext_ni_h(t,u,nocc,nimp,nsites)
#        hop,himp,hcs_imp,hds_imp,hcs_ext,hds_ext,hext=models_embed.hubbard_imp_ext_h(hall_sp,u,nocc,nimp,nsites)

        # g.s embedding basis
        pemb,pcore,pvirt=embed.embed(hlat,nimp,nocc)

        # perturbation operator
        perturb=utils.zeros([nsites,nsites])
        perturb[0,0]=1.

        p_coeffs={}
        p_coeffs[0,0]=1.
        p_coeffs[1,1]=1.
        perturbop=models.ContractedCD(p_coeffs)

        fd=file("ph_siam.out."+str(u),"w")
        for omega in N.arange(0,8,0.1):
    
            ops_dict=response_embed.mb_ph_ops(hlat,perturb,omega+1j*delta,nimp,nocc,pemb,pcore,pvirt)
            configs_dict=response_embed.mb_configs(nsites,nimp,nimp_sp,2*nocc-nimp_sp,0)
            neutral_ops_configs=response_embed.mb_ph_ops_configs(ops_dict,configs_dict)

            # basis is setup, now build matrix representations
            perturb_mat=opbasis_ni.oimp_matrix_form(perturbop,neutral_ops_configs,cocc,vocc)

            # h, neutral configs
            himp_mat=opbasis_ni.oimp_matrix_form(himp,neutral_ops_configs,cocc,vocc)
            himp_ext_mat=opbasis_ni.himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,neutral_ops_configs,cocc,vocc)
            hext_mat=opbasis_ni.hext_matrix_form(hext,neutral_ops_configs,cocc,vocc,e0)
            
            hmat=himp_mat+himp_ext_mat+hext_mat

            unit_mat=opbasis_ni.oimp_matrix_form(models.Unit(),neutral_ops_configs,cocc,vocc)

            # get neutral ground-state energy
            es,cs=la.peigh(hmat,unit_mat)
            e0=es[0]
            psi0=cs[:,0]
            # all matrices setup, solve linear response (complex)
            psi1=la.solve_perturb(omega-1j*delta,unit_mat,hmat,e0,psi0,perturb_mat,project=True)

            ph_gf=N.dot(N.conj(psi1),N.dot(perturb_mat,psi0))

            ph_gf0=_second_order_energy(hlat,nocc,perturb,omega+1j*delta)

            #print omega, ph_gf.imag,ph_gf0.imag
            print omega, ph_gf,2*ph_gf0
            print >>fd, omega-mu, ph_gf.imag,ph_gf0.imag


def sp_greens():
    # main loop for single-particle GF
    utils.dtype=N.complex128

    nimp=1
    nimp_sp=2*nimp
    nocc=20
    nsites=40
    nsites_sp=nsites*2
    ndim=2
    sz=0
    # sites numbered in order as imp+bath ... ext
    sites_imp=range(0,2*nimp_sp)
    sites_ext=range(2*nimp_sp,nsites_sp)
    ncore=2*nocc-2*nimp
    cocc=range(2*nimp_sp,2*nimp_sp+ncore)
    vocc=range(2*nimp_sp+ncore,nsites_sp)

    N.set_printoptions(precision=3)

    t=-1. # hopping
    delta=0.2 # broadening

#    for u in [0.0,4.0,10.0]:
    for u in [1.0]:

        # Single impurity Anderson model
        mu=0.
        hlat,hall,hlat_sp,hall_sp,e0=models_embed.si_anderson_imp_ext_ni_h(t,nocc,nimp,nsites)
        hop,himp,hcs_imp,hds_imp,hcs_ext,hds_ext,hext=models_embed.si_anderson_imp_ext_h(hall_sp,u,nocc,nimp,nsites)

#        single site Hubbard in DMET basis
#        mu=u/2
#        hlat,hall,hlat_sp,hall_sp,e0=models_embed.hubbard_imp_ext_ni_h(t,u,nocc,nimp,nsites)
#        hop,himp,hcs_imp,hds_imp,hcs_ext,hds_ext,hext=models_embed.hubbard_imp_ext_h(hall_sp,u,nocc,nimp,nsites)

        # g.s embedding basis
        pemb,pcore,pvirt=embed.embed(hlat,nimp,nocc)

        # perturbation operator
        perturb=utils.zeros([nsites])
        perturb[0]=1.
        p_coeffs=N.array([1])
        perturb_dop=models.ContractedD(N.array([0]),p_coeffs)
        perturb_cop=models.ContractedC(N.array([0]),p_coeffs)

        fd=file("ph_siam.out."+str(u),"w")
        for omega in N.arange(-2,2,0.1):
    
            ops_dict=response_embed.mb_sp_ops(hlat,perturb,omega+1j*delta,nimp,nocc,pemb,pcore,pvirt)
            configs_dict=response_embed.mb_configs(nsites,nimp,nimp_sp,2*nocc-nimp_sp,0)
            neutral_ops_configs,plus_ops_configs,minus_ops_configs=response_embed.mb_sp_ops_configs(ops_dict,configs_dict)

            # basis is setup, now build matrix representations
            perturb_dop_mat=opbasis_ni.oimp_matrix_bra_ket_form(perturb_dop,minus_ops_configs,neutral_ops_configs,cocc,vocc)
            perturb_cop_mat=opbasis_ni.oimp_matrix_bra_ket_form(perturb_cop,plus_ops_configs,neutral_ops_configs,cocc,vocc)

            # h, N-1 configs
            himp_mat_minus=opbasis_ni.oimp_matrix_form(himp,minus_ops_configs,cocc,vocc)
            himp_ext_mat_minus=opbasis_ni.himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,minus_ops_configs,cocc,vocc)
            hext_mat_minus=opbasis_ni.hext_matrix_form(hext,minus_ops_configs,cocc,vocc,e0)
            hmat_minus=himp_mat_minus+himp_ext_mat_minus+hext_mat_minus
            unit_mat_minus=opbasis_ni.oimp_matrix_form(models.Unit(),minus_ops_configs,cocc,vocc)

            # h, N+1 configs
            himp_mat_plus=opbasis_ni.oimp_matrix_form(himp,plus_ops_configs,cocc,vocc)
            himp_ext_mat_plus=opbasis_ni.himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,plus_ops_configs,cocc,vocc)
            hext_mat_plus=opbasis_ni.hext_matrix_form(hext,plus_ops_configs,cocc,vocc,e0)
            hmat_plus=himp_mat_plus+himp_ext_mat_plus+hext_mat_plus
            unit_mat_plus=opbasis_ni.oimp_matrix_form(models.Unit(),plus_ops_configs,cocc,vocc)

            # h, neutral configs
            himp_mat=opbasis_ni.oimp_matrix_form(himp,neutral_ops_configs,cocc,vocc)
            himp_ext_mat=opbasis_ni.himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,neutral_ops_configs,cocc,vocc)
            hext_mat=opbasis_ni.hext_matrix_form(hext,neutral_ops_configs,cocc,vocc,e0)
            hmat=himp_mat+himp_ext_mat+hext_mat
            unit_mat=opbasis_ni.oimp_matrix_form(models.Unit(),neutral_ops_configs,cocc,vocc)

            # get neutral ground-state energy
            es,cs=la.peigh(hmat,unit_mat)
            e0=es[0]
            psi0=cs[:,0]

            print e0
            print psi0
            # all matrices setup, solve linear response (complex)
            psi1_minus=la.solve_perturb(omega+1j*delta,unit_mat_minus,hmat_minus,e0,psi0,perturb_dop_mat,
                                        project=False,sign_ham=-1.)

            gfminus=N.dot(psi1_minus,N.dot(perturb_dop_mat,psi0))

            psi1_plus=la.solve_perturb(omega+1j*delta,unit_mat_plus,hmat_plus,e0,psi0,perturb_cop_mat,project=False)
        
            gfplus=N.dot(psi1_plus,N.dot(perturb_cop_mat,psi0))

            d_gf,c_gf=_sp_gf(hlat,omega+1j*delta,perturb,nocc)

            print omega, gfplus.imag,c_gf.imag, gfminus.imag,d_gf.imag
            print >>fd, omega-mu, gfplus.imag+gfminus.imag,c_gf.imag+d_gf.imag
