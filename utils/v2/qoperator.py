import numpy as N
import itertools
import utils
zip=itertools.izip

class ProductQOp(object):
    def __init__(self, qops):
        self.qops=qops

    def __call__(self, fermion_config):
        ret=([fermion_config], [1.])
        for qop in reversed(self.qops):
            ret=act(qop, *ret)
        return ret

def act(qop, configs, coeffs):
    ret_configs=[]
    ret_coeffs=[]

    for coeff, config in zip(coeffs, configs):
        new_configs, new_coeffs=qop(config)
        new_coeffs=[nc*coeff for nc in new_coeffs]
        ret_configs+=new_configs
        ret_coeffs+=new_coeffs
    return ret_configs, ret_coeffs
        
# def matrix_form(h, configs, ket_configs=None):
#     if ket_configs==None:
#         ket_configs=configs

#     tuple_configs=[tuple(config) for config in configs]
#     nconfigs=len(configs)
#     configs_indices=dict(zip(tuple_configs, range(nconfigs)))

#     hmat=utils.zeros((nconfigs,len(ket_configs)))
#     for i, config in enumerate(ket_configs):
#         new_configs, hints=h(config)
#         for hint, new_config in zip(hints, new_configs):
#             try: # if new_config is in the set of (bra)configs
#                 hmat[configs_indices[tuple(new_config)],i]+=hint
#             except: # otherwise 0
#                 pass

#     return hmat

def matrix_form(h, bra_configs, ket_configs):
    hmat=utils.zeros([len(bra_configs),len(ket_configs)])

    for i, bra in enumerate(bra_configs):
        for j, ket in enumerate(ket_configs):
            hmat[i,j]=matrix_element(h,bra,ket)
    return hmat

def matrix_element(h, bra_config, ket_config):
    new_configs, hints=h(ket_config)
    return sum([hint for hint, new_config 
                in zip(hints, new_configs) if new_config==bra_config])

def expectation(h, configs, bra, ket):
    tuple_configs=[tuple(config) for config in configs]
    nconfigs=len(configs)
    configs_indices=dict(zip(tuple_configs, range(nconfigs)))

    ret_val=0.
    for i, config in enumerate(configs):
        new_configs, hints=h(config)
        for hint,new_config in zip(hints,new_configs):
            ret_val+=bra[configs_indices[tuple(new_config)]]*hint*ket[i]

    return ret_val
            
