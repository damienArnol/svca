from numpy.core.multiarray import dtype

import numpy as np
import scipy.special

from limix.core.covar import SQExpCov
from limix.core.covar import ZKZCov
from limix.core.covar import FixedCov
from limix.core.covar import SumCov
import limix.core.gp
from limix.core.mean import MeanBase
from limix.utils.preprocess import covar_rescaling_factor_efficient
import sys
import pdb

import pandas as pd


def normalise_phenotype(phenotype):
    phenotype -= phenotype.mean()
    phenotype /= phenotype.std()
    phenotype = np.reshape(phenotype, [len(phenotype), 1])
    return phenotype

def quantile_normalise_phenotype(phenotype):
    # take ranks and scale to uniform
    phenotype_order = pd.Series(phenotype).rank().astype(float)
    phenotype_order = phenotype_order.values
    phenotype_order /= (phenotype_order.max()+1.)

    # transform uniform to gaussian using probit
    phenotype_norm = np.sqrt(2.) * scipy.special.erfinv(2.*phenotype_order-1.)
    phenotype_norm = np.reshape(phenotype_norm, [len(phenotype_norm), 1])

    return phenotype_norm
    # return phenotype_norm


def remove_markers(protein_names, phenotypes):
    markers_to_rm = ['pAKT',
                     'AMPK',
                     'pBAD',
                     'EGFR',
                     'SHP2',
                     'CD45']

    to_rm_ix = np.zeros(len(markers_to_rm))

    for i in range(0, len(markers_to_rm)):
        to_rm_ix[i] = [j for j in range(len(protein_names)) if protein_names[j] == markers_to_rm[i]][0]

    for name in markers_to_rm:
        protein_names.remove(name)

    phenotypes = np.delete(phenotypes, to_rm_ix, 1)

    return protein_names, phenotypes


def build_interaction_matrix(t1, t2, cell_types):
    types_1 = (cell_types == t1) * 1.0
    types_1 = np.reshape(types_1, [len(types_1), 1])

    if t2 == 'all':
        types_2 = np.ones([1, len(cell_types)])
    else:
        types_2 = (cell_types == t2) * 1.0
        types_2 = np.reshape(types_2, [1, len(types_2)])

    return types_1.dot(types_2)


def build_cell_type_kinship(cell_types_rows, cell_types_columns=None):
    if cell_types_columns is None:
        cell_types_columns = cell_types_rows
    cell_types_list = np.unique(cell_types_rows)
    cell_types_list = np.concatenate((np.unique(cell_types_columns), cell_types_list))
    cell_types_list = np.unique(cell_types_list)

    K = None

    for t in cell_types_list:
        cells_selection_columns = np.reshape((cell_types_columns == t), [1, len(cell_types_columns)])
        cells_selection_rows = np.reshape((cell_types_rows == t), [len(cell_types_rows), 1])
        K_tmp = cells_selection_rows.dot(cells_selection_columns)
        if K is None:
            K = K_tmp
        else:
            K += K_tmp

    return K*1.

# TODO clean so that we dont store lists of kernels independently from limix
# TODO usage of int cell type or string should also be more consistent
def build_model(Kinship, phenotype, N_cells, X, cell_types, test_set=None, intrinsic =True, environment = True,
                environmental_cell_types=None, affected_cell_type=None, by_effective_type=True):

    if test_set is not None and test_set.dtype == bool:
        X_training = X[~test_set, :]
        X_test = X[test_set, :]
        mean_training = phenotype[~test_set]
        N_cells_training = sum(~test_set)
        N_cells_test = N_cells - N_cells_training
        cell_types_training = cell_types[~test_set]
        cell_types_test = cell_types[test_set]

    else:
        X_training = X
        X_test = None
        mean_training = phenotype
        N_cells_training = N_cells
        N_cells_test = 0
        cell_types_training = cell_types

    rm_diag = True
    cov = None

    # list cell types
    cell_type_list = np.unique(cell_types)

    # local_noise
    local_noise_cov = SQExpCov(X_training, Xstar=X_test)
    local_noise_cov.setPenalty(mu=50., sigma=50.)

    # noise
    noise_covs = [None for i in range(len(cell_type_list))]
    for t in cell_type_list:
        cells_selection = (cell_types_training == t) * np.eye(N_cells_training)
        if N_cells_test == 0:
            Kcross = None
        # TODO: adapt to multiple cell types
        else:
            # Kcross = np.concatenate((np.zeros([N_cells_test, N_cells_training]), np.eye(N_cells_test)), axis=1)
            Kcross = np.zeros([N_cells_test, N_cells_training])
        noise_covs[t] = FixedCov(cells_selection, Kcross)
        if cov is None:
            cov = SumCov(local_noise_cov, noise_covs[t])
        else:
            cov = SumCov(cov, noise_covs[t])

    # environment effect: for each pair of cell types
    # t1 is the receiving type, t2 is the effector
    if environment:
        if by_effective_type:
            # env_covs = np.array([len(cell_type_list), len(cell_type_list)])
            env_covs = [[None for i in range(len(cell_type_list))] for j in range(len(cell_type_list))]
        else:
            env_covs = [None for i in range(len(cell_type_list))]
        # env_covs = [tmp] * len(cell_type_list)
        for t1 in cell_type_list:
            if affected_cell_type is not None and affected_cell_type != t1:
                continue
            if by_effective_type:
                for t2 in cell_type_list:
                    # select only the environmental cell type if not all
                    if environmental_cell_types is not None and environmental_cell_types != t2:
                        continue
                    interaction_matrix = build_interaction_matrix(t1, t2, cell_types)
                    tmp = ZKZCov(X, Kinship, rm_diag, interaction_matrix, test_set)
                    env_covs[t1][t2] = tmp
                    env_covs[t1][t2].setPenalty(mu=200., sigma=50.)
                    cov = SumCov(cov, env_covs[t1][t2])
            else:
                interaction_matrix = build_interaction_matrix(t1, 'all', cell_types)
                tmp = ZKZCov(X, Kinship, rm_diag, interaction_matrix, test_set)
                env_covs[t1] = tmp
                env_covs[t1].setPenalty(mu=200., sigma=50.)
                cov = SumCov(cov, env_covs[t1])
    else:
        env_covs = None

    if intrinsic:
        K = build_cell_type_kinship(cell_types_training)
        if N_cells_test != 0:
            Kcross = build_cell_type_kinship(cell_types_test, cell_types_training)
        intrinsic_cov = FixedCov(K, Kcross)
        cov = SumCov(cov, intrinsic_cov)
    else:
        intrinsic_cov = None

    # mean term
    mean = MeanBase(mean_training)

    # define GP
    gp = limix.core.gp.GP(covar=cov, mean=mean)

    print 'GP created '

    return gp, noise_covs, local_noise_cov, env_covs, intrinsic_cov

def add_env_to_null(gp_init, cell_types, phenotype, X, env_type_ix,
                    affected_type_ix):
    cov = gp_init.covar

    Kinship = build_cell_type_kinship(cell_types)
    cell_type_list = np.unique(cell_types)
    rm_diag=True

    env_covs = [[None for i in range(len(cell_type_list))]for j in range(len(cell_type_list))]
    for t1 in cell_type_list:
        if affected_type_ix is not None and affected_type_ix != t1:
            continue
        for t2 in cell_type_list:
            # select only the environmental cell type if not all
            if env_type_ix is not None and env_type_ix != t2:
                continue
            interaction_matrix = build_interaction_matrix(t1, t2, cell_types)
            env_covs[t1][t2] = ZKZCov(X, Kinship, rm_diag, interaction_matrix)
            env_covs[t1][t2].setPenalty(mu=200, sigma=50)
            cov = SumCov(cov, env_covs[t1][t2])

    gp_init.covar = cov
    gp_init.optimize()
    return gp_init

def compute_r2(truth, prediction):
    ss_tot = ((truth - truth.mean()) ** 2.0).sum()
    ss_res = ((truth - prediction) ** 2.0).sum()
    r2 = 1.0 - ss_res / ss_tot
    return r2


def read_data(expression_file, position_file, protein_index=None):
    with open(expression_file, 'r') as f:
        prot_tmp = f.readline()

    protein_names = prot_tmp.split(' ')
    protein_names[-1] = protein_names[-1][0:-1]  # removing the newline sign at the end of the last protein
    phenotypes = np.loadtxt(expression_file, delimiter=' ', skiprows=1)

    # TODO sort out
    #protein_names, phenotypes = remove_markers(protein_names, phenotypes)

    protein_names = np.reshape(protein_names, [len(protein_names), 1])

    X = np.genfromtxt(position_file, delimiter=',')

    if X.shape[0] != phenotypes.shape[0]:
        raise Exception('cell number inconsistent between position and expression levels ')

    if protein_index is not None:
        protein_names = protein_names[protein_index, :]
        phenotypes = phenotypes[:, protein_index]

    return protein_names, phenotypes, X


# cell_types_names is a list of unique names
def read_types(cell_types_file, with_names = False):
    cell_types_list = np.loadtxt(cell_types_file, dtype=str)
    cell_types = np.zeros(len(cell_types_list), dtype=int)
    cell_types_names = np.unique(cell_types_list)
    n_types = len(cell_types_names)

    for t in range(0, n_types):
        cell_types[cell_types_list == cell_types_names[t]] = t

    if with_names:
        return cell_types, cell_types_names, n_types, cell_types_list
    return cell_types, cell_types_names, n_types

# TODO decide whether to use parameters directly or Gower normalisation (not simple when considering the targeter types separately)
def write_results(trained_model, output_dir, file_prefix, n_types, cell_types, cell_types_names, by_effective_type=True):
    # extracting covariance terms from the trained model
    ####################################################################
    intrinsic_cov = trained_model['intrinsic_cov']
    noise_covs = trained_model['noise_covs']
    local_noise_cov = trained_model['local_noise_cov']
    env_covs = trained_model['env_covs']

    # putting parameters together
    ####################################################################
    if by_effective_type:
        parameters = np.zeros([1, int(2 * n_types**2.0 + n_types + 2 + 1)])
    else:
        parameters = np.zeros([1, int(2 * n_types + n_types + 2 + 1)])

    try:
        parameters[0,0] = 1./covar_rescaling_factor_efficient(intrinsic_cov.K())
    except:
        parameters[0,0] = np.nan
    try:
        parameters[0,1] = 1./covar_rescaling_factor_efficient(local_noise_cov.K())
    except:
        parameters[0,1] = np.nan
    parameters[0,2] = local_noise_cov.length
    count = 3

    for type1 in np.unique(cell_types):
        cell_filter = (cell_types == type1)
        assert np.all(noise_covs[type1].K()[:, ~cell_filter][~cell_filter,:] == 0), 'problem cell filter'
        K_tmp = noise_covs[type1].K()[:, cell_filter][cell_filter,:]
        try:
            parameters[0,count] = 1./covar_rescaling_factor_efficient(K_tmp)
        except:
            parameters[0,count] = np.nan
        count += 1

    for type1 in np.unique(cell_types):
        cell_filter = (cell_types == type1)
        if by_effective_type:
            for type2 in np.unique(cell_types):
                K_tmp = env_covs[type1][type2].K()[:, cell_filter][cell_filter,:]
                assert np.all(env_covs[type1][type2].K()[:, ~cell_filter][~cell_filter,:]  ==0), 'problem cell filter'
                try:
                    parameters[0,count] = 1./covar_rescaling_factor_efficient(K_tmp)
                except:
                    parameters[0,count] = np.nan
                count += 1
        else:
            K_tmp = env_covs[type1].K()[:, cell_filter][cell_filter,:]
            assert np.all(env_covs[type1].K()[:, ~cell_filter][~cell_filter,:]  ==0), 'problem cell filter'
            try:
                parameters[0,count] = 1./covar_rescaling_factor_efficient(K_tmp)
            except:
                parameters[0,count] = np.nan
            count += 1

    for type1 in np.unique(cell_types):
        if by_effective_type:
            for type2 in np.unique(cell_types):
                parameters[0,count] = env_covs[type1][type2].length
                count += 1
        else:
            parameters[0,count] = env_covs[type1].length
            count += 1

    # Putting header together
    ####################################################################
    if by_effective_type:
        result_header_array = [None for i in range(int(2 * n_types**2.0 + n_types + 2 + 1))]
    else:
        result_header_array = [None for i in range(int(2 * n_types + n_types + 2 + 1))]
    result_header_array[0] = 'inrtinsic'
    result_header_array[1] = 'local_noise_scale'
    result_header_array[2] = 'local_noise_length'
    count = 3

    for type in np.unique(cell_types):
            result_header_array[count] = 'noise_'+cell_types_names[type]
            count += 1

    for type1 in np.unique(cell_types):
        if by_effective_type:
            for type2 in np.unique(cell_types):
                result_header_array[count] = 'env_'+cell_types_names[type1]+'_'+cell_types_names[type2]
                count += 1
        else:
            result_header_array[count] = 'env_'+cell_types_names[type1]
            count += 1

    for type1 in np.unique(cell_types):
        if by_effective_type:
            for type2 in np.unique(cell_types):
                result_header_array[count] = 'env_length_'+cell_types_names[type1]+'_'+cell_types_names[type2]
                count += 1
        else:
            result_header_array[count] = 'env_length_'+cell_types_names[type1]
            count += 1

    # TODO concatenate list into single string
    result_header = ' '.join(result_header_array)

    # Write file
    ####################################################################
    output_file = output_dir + '/' + file_prefix
    with open(output_file, 'w') as f:
        np.savetxt(f,
                   parameters,
                   delimiter=' ',
                   header=result_header,
                   fmt='%s',
                   comments='')
    # log_lik_file = output_file + '_loglik'
    # with open(log_lik_file, 'w') as f:
    #     np.savetxt(f, log_lik)

def subset_image(X, phenotype, cell_types, types_by_name, size = 500, order=False):
    sel = np.logical_and(X[:,0] < size, X[:,1] < size)

    X = X[sel, :]
    phenotype = phenotype[sel]
    cell_types = cell_types[sel]
    types_by_name = types_by_name[sel]

    if order:
        c_order = np.argsort(cell_types)

        X = X[c_order,:]
        phenotype = phenotype[c_order]
        cell_types = cell_types[c_order]
        types_by_name = types_by_name[c_order]

    return X, phenotype, cell_types, types_by_name
