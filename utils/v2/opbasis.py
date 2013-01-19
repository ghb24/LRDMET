import collections
import numpy as N
import qoperator
import utils

"""
Functions to create operator representations in
an "opbasis". An opbasis is a set of configs multiplied by a set of ops.
Thus matrix elements are of the form <basis opi| qop | opj basis>.
These functions are all slow but general. Specialized hard-coded forms
for various operators are in opbasis_ni.
"""
def matrix_bra_ket_form(qop,ops_bra_configs,ops_ket_configs):
    """
    < bra opi | qop | opj ket >
    where <bra opi| and |opj ket> can be different sets
    """
    bra_nconfigs=[len(configs) for (ops, configs) in ops_bra_configs]
    ket_nconfigs=[len(configs) for (ops, configs) in ops_ket_configs]
    full_mat=utils.zeros([sum(bra_nconfigs),sum(ket_nconfigs)])

    for i, (opi, bra_configsi) in enumerate(ops_bra_configs):
        iptr=sum(bra_nconfigs[:i])
        for j, (opj, ket_configsj) in enumerate(ops_ket_configs):
            jptr=sum(ket_nconfigs[:j])
            prod_op=qoperator.ProductQOp([opi.t(),qop,opj])
            mat=qoperator.matrix_form(prod_op,bra_configsi,ket_configsj)
            full_mat[iptr:iptr+bra_nconfigs[i],jptr:jptr+ket_nconfigs[j]]=mat
    return full_mat

def matrix_form(qop,ops_configs):
    """
    < bra opi | qop | opj ket >
    where <bra opi| and |opj ket> are the same set
    """
    nconfigs=[len(configs) for (ops, configs) in ops_configs]
    full_mat=utils.zeros([sum(nconfigs),sum(nconfigs)])

    for i, (opi, bra_configsi) in enumerate(ops_configs):
        iptr=sum(nconfigs[:i])
        for j, (opj, ket_configsj) in enumerate(ops_configs):
            jptr=sum(nconfigs[:j])
            prod_op=qoperator.ProductQOp([opi.t(),qop,opj])
            mat=qoperator.matrix_form(prod_op,bra_configsi,ket_configsj)
            full_mat[iptr:iptr+nconfigs[i],jptr:jptr+nconfigs[j]]=mat
    return full_mat

def matrix_oneside_form(qop,ops_configs):
    """
    alternative algorithm to matrix_form
    """
    nconfigs=[len(configs) for (ops, configs) in ops_configs]
    full_mat=utils.zeros([sum(nconfigs),sum(nconfigs)])

    for i, (opi, bra_configsi) in enumerate(ops_configs):
        iptr=sum(nconfigs[:i])
        for j, (opj, ket_configsj) in enumerate(ops_configs):
            jptr=sum(nconfigs[:j])
            mat=matrix_opop_element(qop, opi, bra_configsi, opj, ket_configsj)
            full_mat[iptr:iptr+nconfigs[i],jptr:jptr+nconfigs[j]]=mat

    return full_mat

def matrix_opop_element(h, bra_op, bra_configs, ket_op, ket_configs):
    hmat=utils.zeros([len(bra_configs),len(ket_configs)])
    hket_op=qoperator.ProductQOp([h,ket_op])
    for i, bra_config in enumerate(bra_configs):
        new_bra_configs, bints=bra_op(bra_config)
        bra_dict=collections.defaultdict(float)
        for nbc, bi in zip(new_bra_configs,bints):
            bra_dict[tuple(nbc)]+=bi
        
        for j, ket_config in enumerate(ket_configs):
            new_ket_configs, kints=hket_op(ket_config)

            for nkc, kj in zip(new_ket_configs,kints):
                try:
                    hmat[i,j]+=bra_dict[tuple(nkc)]*kj
                except:
                    pass

    return hmat
