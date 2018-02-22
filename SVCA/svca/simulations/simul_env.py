import numpy as np
from svca.models.model1 import Model1
from svca.simulations.from_real import FromRealSimulation
from svca.models.io import *
from svca.util_functions import utils
import sys
from limix.utils.preprocess import covar_rescaling_factor_efficient
from copy import deepcopy


def run(data_dir, protein_index, output_dir, interactions_size,
        normalisation='standard', permute=False):
    # reading all data
    ####################################################################
    expression_file = data_dir + '/expressions.txt'
    position_file = data_dir+'/positions.txt'
    protein_names, phenotypes, X = utils.read_data(expression_file,
                                                   position_file)

    protein_name = protein_names[protein_index, :]
    phenotype = phenotypes[:, protein_index]
    sel = range(phenotypes.shape[1])
    sel.remove(protein_index)
    kin_from = phenotypes[:, sel]

    # N_samples = X.shape[0]

    boot_ix = deepcopy(interactions_size)
    interactions_size = float(int(interactions_size)%10)/10.
    down_sampling = 1 - float(int(interactions_size)/10)/10.

    # down sampling
    n_sel = down_sampling * X.shape[0]
    sel = np.sort(np.random.choice(range(X.shape[0]), n_sel, replace=False))
    X = X[sel,:]
    phenotype = phenotype[sel]
    kin_from= kin_from[sel,:]
    N_samples = X.shape[0]
    # TODO select X, select phenotype, kin_from, N_samples


    # permuting cells
    if permute:
        perm = np.random.permutation(X.shape[0])
        X = X[perm, :]

    # do null simulation
    ####################################################################
    sim = FromRealSimulation(X, phenotype, kin_from)
    Y_sim = sim.simulate(interactions_size=interactions_size)

    # run model on simulated data
    ####################################################################
    # all but interactions
    ####################################################################
    cterms = ['intrinsic', 'environmental']
    model = Model1(Y_sim, X, norm=normalisation, oos_predictions=0., cov_terms=cterms, kin_from=kin_from)
    model.reset_params()
    model.train_gp(grid_size=10)

    file_prefix = protein_name[0] + '_' + str(boot_ix) + '_environmental'
    write_variance_explained(model, output_dir, file_prefix)
    write_LL(model, output_dir, file_prefix)

    ####################################################################
    # adding interactions
    ####################################################################
    model.add_cov(['interactions'])
    model.reset_params()
    model.train_gp(grid_size=10)

    file_prefix = protein_name[0] + '_' + str(boot_ix) + '_interactions'
    write_variance_explained(model, output_dir, file_prefix)
    write_LL(model, output_dir, file_prefix)


if __name__ == '__main__':
    data_dir = sys.argv[1]
    output_dir = sys.argv[2]
    protein_index = int(sys.argv[3])
    bootstrap_index = sys.argv[4]
    interactions_size = bootstrap_index
    normalisation = sys.argv[5]


    # data_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/IMC_reprocessed_median/Cy1x7/'
    # # protein_index = 19
    # protein_index = 3 # 5
    # interactions_size = 3
    # output_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/tests/test_new_init/'
    # normalisation='quantile'
    # perm = False

    run(data_dir, protein_index, output_dir, interactions_size, normalisation)
