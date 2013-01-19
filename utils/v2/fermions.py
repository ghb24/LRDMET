"""
Functions to deal with fermion configurations

|config> is in occupation number representation
list of 0, 1s etc. of length 2*number of physical sites

By convention, EVEN sites are alpha sites
               ODD sites are beta sites
"""

import numpy as N
import random
import choose

def destroy(config, pos):
    if config==[]:
        return 0, []

    if config[pos]:
        new_config=config[:]         
        new_config[pos]=0
        n=sum(config[:pos])
        sign=-1 if n & 1 else 1
        return sign, new_config
    else:
        return 0, []

def create(config, pos):
    if config==[]:
        return 0, []

    if not config[pos]:
        new_config=config[:] 
        new_config[pos]=1
        n=sum(config[:pos])
        sign=-1 if n & 1 else 1
        return sign, new_config
    else:
        return 0, []

def excite(config, first, second):
    """
    apply a+(first) a(second)|config>
    return phase, new_config
    phase=0, new_config=[] if annihilated
    """
    if config==[]:
        return 0, []

    if first==second and config[second]:
        return 1, config[:]
    elif config[first]:
        # cannot excite to since first is occupied 
        return 0, []
    elif not config[second]:
        # cannot excite from since second spot is occupied
        return 0, []
    else:    
        new_config=config[:]
        new_config[second], new_config[first]=0, 1
        return phase(config, first, second), new_config

def phase(config, first, second):
    """
    phase associated with a+(first) a(second)|config>
    """
    lower, upper=min(first, second), max(first, second)
    nswaps=sum(config[lower+1:upper])
    return -1 if nswaps & 1 else 1

def to_string(occs, nsites):
    occ_string=N.zeros([nsites],N.int)
    occ_string[occs]=1
    return list(occ_string)

def to_occs(config):
    """
    occ representation -> labelled representation
    'string' is a seq of 1, 0's.
    returns index of all '1's in string
    """
    return [index for index, el in enumerate(config) if el]


def spinless_random_conserving_config(nsites, totaln):
    assert 0<=totaln<=nsites
    occ_indices=set()
    while len(occ_indices)<totaln:
        occ_indices.add(int(N.random.randint(nsites)))
    return to_string(list(occ_indices),nsites)


def configs(nsites, ndim):
    """
    int -> list(int)
    Generate list of all possible fermion configs
    associated with a given number of sites
    """
    return [list(config) for config in 
            N.ndindex(tuple([ndim]*nsites))]

def conserving_configs(nsites, ndim, totaln, totalsz):
    """
    Generate list of all possible fermion
    configs with given sz and given N
    """
    nalpha=(totaln+totalsz)/2
    nbeta=totaln-nalpha
    assert nalpha <= nsites/2 and nbeta <=nsites/2
    configs=[]
    for alpha_occ in choose.choose(range(nsites/2), nalpha):
        for beta_occ in choose.choose(range(nsites/2), nbeta):
            configs.append(assemble(to_string(alpha_occ, nsites/2),
                                    to_string(beta_occ, nsites/2)))
    return configs

def spinless_conserving_configs(nsites, totaln):
    """
    Generate list of all possible fermion
    configs with given N
    """
    return [to_string(choice, nsites) for choice
            in choose.choose(range(nsites), totaln)]

def flip(spin):
    return (spin+1) & 1

def alpha_string(config):
    return config[0:len(config):2]

def beta_string(config):
    return config[1:len(config):2]

def assemble(alpha_string, beta_string):
    config=[None]*2*len(alpha_string)
    config[0:len(config):2]=alpha_string
    config[1:len(config):2]=beta_string
    return config
