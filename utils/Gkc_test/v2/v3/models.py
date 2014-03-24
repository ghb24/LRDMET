import fermions
import qoperator

import numpy as N
import scipy.linalg

"""
operator \sum_i c_i a+_i 
"""
class ContractedC(object):
    def __init__(self, sites, site_coeffs):
        self.sites=sites[:]
        self.site_coeffs=site_coeffs[:]
        self.fermion=True

    def __call__(self, fermion_config):
        new_configs=[]
        new_vals=[]
        for site, site_coeff in zip(self.sites, self.site_coeffs):
            sign, new_config=fermions.create(fermion_config, site)
        
            if sign:
                new_configs.append(new_config)
                new_vals.append(site_coeff*sign)
        return new_configs, new_vals
    
    def coeffs(self, i):
        try:
            ix=self.sites.index(i)
            return self.site_coeffs[ix]
        except:
            print i, self.sites
            assert False
                               
    def __str__(self):
        return "C:"+str(self.sites)+str(self.site_coeffs)

    def t(self):
        return ContractedD(self.sites, N.conj(self.site_coeffs))

"""
operator \sum_i c_i a_i 
"""
class ContractedD(object):
    def __init__(self, sites, site_coeffs):
        self.sites=sites[:]
        self.site_coeffs=site_coeffs[:]
        self.fermion=True

    def __call__(self, fermion_config):
        new_configs=[]
        new_vals=[]
        for site, site_coeff in zip(self.sites, self.site_coeffs):
            sign, new_config=fermions.destroy(fermion_config, site)
        
            if sign:
                new_configs.append(new_config)
                new_vals.append(site_coeff*sign)
        return new_configs, new_vals

    def coeffs(self, i):
        try:
            ix=self.sites.index(i)
            return self.site_coeffs[ix]
        except:
            print i, self.sites
            assert False

    def __str__(self):
        return "D:"+str(self.sites)+str(self.site_coeffs)

    def t(self):
        return ContractedC(self.sites, N.conj(self.site_coeffs))

"""
operator \sum_i c_ij a+i a_j
"""
class ContractedCD(object):
    def __init__(self, coeffs, thresh=1.e-12):
        self.coeffs_dict=coeffs
        self.thresh=thresh
        self.fermion=False

    def __call__(self, fermion_config):
        new_configs=[]
        new_vals=[]
#        for index in N.ndindex(self.coeffs.shape):
#            if abs(self.coeffs[index])>self.thresh:
        for index in self.coeffs_dict:
            sign, new_config=fermions.excite(fermion_config, index[0], index[1])

            if sign:
                new_configs.append(new_config)
                new_vals.append(self.coeffs_dict[index[0],index[1]]*sign)

        return new_configs, new_vals

    def coeffs(self, i, j):
        try:
            return self.coeffs_dict[i,j]
        except KeyError:
            return 0.

    def __str__(self):
        return "CD:"+str(self.coeffs_dict)

    def t(self):
#        new_coeffs=self.coeffs.T
        new_coeffs={}
        for key in self.coeffs_dict:
            new_coeffs[tuple(reversed(key))]=N.conj(self.coeffs_dict[key])
        return ContractedCD(new_coeffs)

class Unit(object):

    def __init__(self):
        self.fermion=False

    def __call__(self, fermion_config):
        return [fermion_config], [1.]

    def __str__(self):
        return "Unit"

    def t(self):
        return Unit()

"""
\sum_{ij} t_{ij} a+i a_j + \sum_i U_i nia nib
"""
class GeneralHubbard(object):

    def __init__(self, t_dict, u_dict):
        self.t_dict=t_dict.copy()
        self.u_dict=u_dict.copy()
        self.fermion=False

    def __call__(self, fermion_config):
        new_configs=[]
        new_vals=[]
        for index in self.t_dict.keys():
            sign, new_config=fermions.excite(fermion_config, index[0], index[1])
        
            if sign:
                new_configs.append(new_config)
                new_vals.append(self.t_dict[index[0],index[1]]*sign)

        for pair in self.u_dict.keys():
            if fermion_config[pair[0]] and fermion_config[pair[1]]:
                new_configs.append(fermion_config[:])
                new_vals.append(self.u_dict[pair[0],pair[1]])

        return new_configs, new_vals

    def __str__(self):
        return "CD:"+str(self.t_dict)+"U"+str(self.u_dict)

    def t(self):
        return GeneralHubbard(self.t_dict, self.u_dict)

