import scipy.linalg
import numpy as N
import models
import qoperator
import utils

def eval_hd_cv(hd_ext,opj,cocc,vocc):
    val=0.
    for v in vocc:
        val+=hd_ext.coeffs(v)*opj.coeffs(v)
    return val

def eval_codv_hd_cv(opit,hd_ext,opj,cocc,vocc):
    val=0.
    for c in cocc:
        for v in vocc:
          val-=opit.coeffs(c,v)*hd_ext.coeffs(c)*opj.coeffs(v)
    return val

def eval_co_hd(opit,hd_ext,cocc,vocc):
    val=0.
    for c in cocc:
        val+=opit.coeffs(c)*hd_ext.coeffs(c)
    return val

def eval_co_hd_cvdo(opit,hd_ext,opj,cocc,vocc):
    val=0.
    for c in cocc:
        for v in vocc:
            val+=opit.coeffs(c)*hd_ext.coeffs(v)*opj.coeffs(v,c)
    return val

def eval_hc_do(hc_ext,opj,cocc,vocc):
    val=0.
    for c in cocc:
        val+=hc_ext.coeffs(c)*opj.coeffs(c)
    return val
    
def eval_codv_hc_do(opit,hc_ext,opj,cocc,vocc):
    val=0.
    for c in cocc:
        for v in vocc:
          val+=opit.coeffs(c,v)*hc_ext.coeffs(v)*opj.coeffs(c)
    return val

def eval_dv_hc(opit,hc_ext,cocc,vocc):
    val=0.
    for v in vocc:
        val+=opit.coeffs(v)*hc_ext.coeffs(v)
    return val
    
def eval_dv_hc_cvdo(opit,hc_ext,opj,cocc,vocc):
    val=0.
    for c in cocc:
        for v in vocc:
            val-=opit.coeffs(v)*hc_ext.coeffs(c)*opj.coeffs(v,c)
    return val

def eval_h(hext,cocc,vocc,e0):
    return e0

def eval_dv_h_cv(opit,hext,opj,cocc,vocc,e0):
    val=0.
    for v1 in vocc:
        for v2 in vocc:
            val+=(hext.coeffs(v1,v2))*opit.coeffs(v1)*opj.coeffs(v2)
    for v1 in vocc:
        val+=e0*opit.coeffs(v1)*opj.coeffs(v1)
    return val

def eval_co_h_do(opit,hext,opj,cocc,vocc,e0):
    val=0.
    for c1 in cocc:
        for c2 in cocc:
            val-=hext.coeffs(c1,c2)*opit.coeffs(c1)*opj.coeffs(c2)
    for c1 in cocc:
        val+=e0*opit.coeffs(c1)*opj.coeffs(c1)
    return val

def eval_codv_h_cvdo(opit,hext,opj,cocc,vocc,e0):
    val=0.
    # off-diagonal terms
    for c1 in cocc:
        for v1 in vocc:
            for v2 in vocc:
                if v1!=v2:
                    val+=hext.coeffs(v1,v2)*opit.coeffs(c1,v1)*opj.coeffs(v2,c1)
    for c1 in cocc:
        for c2 in cocc:
            for v1 in vocc:
                if c1!=c2:
                    val+=(-hext.coeffs(c1,c2))*opit.coeffs(c1,v1)*opj.coeffs(v1,c2)
    # diagonal terms
    for c1 in cocc:
        for v1 in vocc:
            val+=(hext.coeffs(v1,v1)-hext.coeffs(c1,c1)+e0)*opit.coeffs(c1,v1)*opj.coeffs(v1,c1)
    return val

def norm_d(opi,opj,cocc,vocc):
    norm=0.
    for c in cocc:
        norm+=opi.coeffs(c)*opj.coeffs(c)
    return norm

def norm_c(opi,opj,cocc,vocc):
    norm=0.
    for v in vocc:
        norm+=opi.coeffs(v)*opj.coeffs(v)
    return norm

def norm_cd(opi,opj,cocc,vocc):
    norm=0.
    for c in cocc:
        for v in vocc:
            norm+=opi.coeffs(c,v)*opj.coeffs(v,c)
    return norm

def oimp_matrix_bra_ket_form(oimp,ops_bra_configs,ops_ket_configs,cocc,vocc):
    # imp is an operator that acts only on the imp and bath orbitals.
    # (e.g. himp)
    bra_start={}
    bra_end={}
    ket_start={}
    ket_end={}
    bra_nconfigs=[len(bra_configs) for (ops, bra_configs) in ops_bra_configs]
    ket_nconfigs=[len(ket_configs) for (ops, ket_configs) in ops_ket_configs]
    bra_ops=[opi for (opi, bra_configs) in ops_bra_configs]
    ket_ops=[opi for (opi, ket_configs) in ops_ket_configs]
    for i, opi in enumerate(bra_ops):
        bra_start[opi]=sum(bra_nconfigs[:i])
        bra_end[opi]=bra_start[opi]+bra_nconfigs[i]

    for i, opi in enumerate(ket_ops):
        ket_start[opi]=sum(ket_nconfigs[:i])
        ket_end[opi]=ket_start[opi]+ket_nconfigs[i]

    full_mat=utils.zeros([sum(bra_nconfigs),sum(ket_nconfigs)])
    
    for i, (opi, bra_configsi) in enumerate(ops_bra_configs):
        opit=opi.t()
        
        for j, (opj, ket_configsj) in enumerate(ops_ket_configs):
            # <imp1.core|opi oimp opj|core.imp2>:
            # parity is opi.parity*oimp.parity*core.parity
            # but assume core is even
            mat=qoperator.matrix_form(oimp,bra_configsi,ket_configsj)
            parity=1.
            if oimp.fermion and opi.fermion:
                parity=-1.

            if isinstance(opit,models.Unit) and isinstance(opj,models.Unit):
                norm=1.
                full_mat[bra_start[opi]:bra_end[opi],ket_start[opj]:ket_end[opj]]=norm*mat*parity
            elif isinstance(opi,models.ContractedC) and isinstance(opj,models.ContractedC):
                norm=norm_c(opit,opj,cocc,vocc)
                full_mat[bra_start[opi]:bra_end[opi],ket_start[opj]:ket_end[opj]]=norm*mat*parity
            elif isinstance(opi,models.ContractedD) and isinstance(opj,models.ContractedD):
                norm=norm_d(opit,opj,cocc,vocc)
                full_mat[bra_start[opi]:bra_end[opi],ket_start[opj]:ket_end[opj]]=norm*mat*parity
            elif isinstance(opi,models.ContractedCD) and isinstance(opj,models.ContractedCD):
                norm=norm_cd(opit,opj,cocc,vocc)
                full_mat[bra_start[opi]:bra_end[opi],ket_start[opj]:ket_end[opj]]=norm*mat*parity

    return full_mat

def oimp_matrix_form(oimp,ops_configs,cocc,vocc):
    # matrix elements of 
    # sum < | oimp | > where oimp acts only on the imp+bath space
    # cocc, vocc: lists of core, virtual labels
    assert not oimp.fermion # operator should have an even number of c/d operators
    start={}
    end={}
    nconfigs=[len(configs) for (ops, configs) in ops_configs]
    ops=[opi for (opi, configi) in ops_configs]
    for i, opi in enumerate(ops):
        start[opi]=sum(nconfigs[:i])
        end[opi]=start[opi]+nconfigs[i]

    full_mat=utils.zeros([sum(nconfigs),sum(nconfigs)])
    
    for i, (opi, configsi) in enumerate(ops_configs):
        opit=opi.t()
        mat=qoperator.matrix_form(oimp,configsi,configsi)
        for j, (opj, configj) in enumerate(ops_configs):
            if isinstance(opit,models.Unit) and isinstance(opj,models.Unit):
                norm=1.
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]=norm*mat
            elif isinstance(opi,models.ContractedC) and isinstance(opj,models.ContractedC):
                norm=norm_c(opit,opj,cocc,vocc)
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]=norm*mat
            elif isinstance(opi,models.ContractedD) and isinstance(opj,models.ContractedD):
                norm=norm_d(opit,opj,cocc,vocc)
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]=norm*mat
            elif isinstance(opi,models.ContractedCD) and isinstance(opj,models.ContractedCD):
                norm=norm_cd(opit,opj,cocc,vocc)
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]=norm*mat


    return full_mat

def himp_ext_matrix_form(hcs_imp,hds_ext,hds_imp,hcs_ext,ops_configs,cocc,vocc):
    # matrix elements of cross terms 
    # sum < |  hcs_imp otimes hds_ext | > + < | hds_imp otimes hcs_ext | >
    assert not len(cocc) & 1, "parities assume even #els in core"

    start={}
    end={}
    nconfigs=[len(configs) for (ops, configs) in ops_configs]
    ops=[opi for (opi, configi) in ops_configs]
    for i, opi in enumerate(ops):
        start[opi]=sum(nconfigs[:i])
        end[opi]=start[opi]+nconfigs[i]

    full_mat=utils.zeros([sum(nconfigs),sum(nconfigs)])

    for i, (opi, configi) in enumerate(ops_configs):
        opit=opi.t()
        for j, (opj, configj) in enumerate(ops_configs):
            for hc_imp, hd_ext in zip(hcs_imp,hds_ext):
                # <a1 b1.core|cimp dext|b2.core a2>: 
                # parity =-1 if fermion(b1) and core is even
                if isinstance(opit,models.Unit) and isinstance(opj,models.ContractedC):
                    mat_imp=qoperator.matrix_form(hc_imp,configi,configj)                    
                    val=eval_hd_cv(hd_ext,opj,cocc,vocc)
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

                elif isinstance(opit,models.ContractedCD) and isinstance(opj,models.ContractedC):
                    mat_imp=qoperator.matrix_form(hc_imp,configi,configj)                    
                    val=eval_codv_hd_cv(opit,hd_ext,opj,cocc,vocc)
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

                # parity is -1, since b1 is fermionic
                elif isinstance(opit,models.ContractedC) and isinstance(opj,models.Unit):
                    mat_imp=qoperator.matrix_form(hc_imp,configi,configj)                    
                    val=-eval_co_hd(opit,hd_ext,cocc,vocc) # parity=-1
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

                elif isinstance(opit,models.ContractedC) and isinstance(opj,models.ContractedCD):
                    mat_imp=qoperator.matrix_form(hc_imp,configi,configj)                    
                    val=-eval_co_hd_cvdo(opit,hd_ext,opj,cocc,vocc) #parity=-1
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

            for hd_imp, hc_ext in zip(hds_imp,hcs_ext):
                # <a1 b1.core|cext dimp|b2.core a2>: 
                # parity =-1 if fermion(b2) and core is even
                if isinstance(opit,models.Unit) and isinstance(opj,models.ContractedD):
                    mat_imp=qoperator.matrix_form(hd_imp,configi,configj)                    
                    val=-eval_hc_do(hc_ext,opj,cocc,vocc) # parity=-1
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

                elif isinstance(opit,models.ContractedCD) and isinstance(opj,models.ContractedD):
                    mat_imp=qoperator.matrix_form(hd_imp,configi,configj)                    
                    val=-eval_codv_hc_do(opit,hc_ext,opj,cocc,vocc) # parity=-1
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

                elif isinstance(opit,models.ContractedD) and isinstance(opj,models.Unit):
                    mat_imp=qoperator.matrix_form(hd_imp,configi,configj)                    
                    val=eval_dv_hc(opit,hc_ext,cocc,vocc) 
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

                elif isinstance(opit,models.ContractedD) and isinstance(opj,models.ContractedCD):
                    mat_imp=qoperator.matrix_form(hd_imp,configi,configj)                    
                    val=eval_dv_hc_cvdo(opit,hc_ext,opj,cocc,vocc) 
                    full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val
                
    return full_mat

def hext_matrix_form(hext,ops_configs,cocc,vocc,e0):
    # matrix elements of 
    # sum < | oext | > where oimp acts only on the imp+bath space
    # cocc, vocc: lists of core, virtual labels
    # e0: core energy
    start={}
    end={}
    nconfigs=[len(configs) for (ops, configs) in ops_configs]
    ops=[opi for (opi, configi) in ops_configs]
    for i, opi in enumerate(ops):
        start[opi]=sum(nconfigs[:i])
        end[opi]=start[opi]+nconfigs[i]

    full_mat=utils.zeros([sum(nconfigs),sum(nconfigs)])

    # hext
    for i, (opi, configi) in enumerate(ops_configs):
        opit=opi.t()
        mat_imp=N.eye(len(configi))
        for j, (opj, configj) in enumerate(ops_configs):
            if isinstance(opit,models.Unit) and isinstance(opj,models.Unit):
                val=eval_h(hext,cocc,vocc,e0)
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

            elif isinstance(opit,models.ContractedD) and isinstance(opj,models.ContractedC):
                val=eval_dv_h_cv(opit,hext,opj,cocc,vocc,e0)
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

            elif isinstance(opit,models.ContractedC) and isinstance(opj,models.ContractedD):
                val=eval_co_h_do(opit,hext,opj,cocc,vocc,e0)
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

            elif isinstance(opit,models.ContractedCD) and isinstance(opj,models.ContractedCD):
                val=eval_codv_h_cvdo(opit,hext,opj,cocc,vocc,e0)
                full_mat[start[opi]:end[opi],start[opj]:end[opj]]+=mat_imp*val

    return full_mat
