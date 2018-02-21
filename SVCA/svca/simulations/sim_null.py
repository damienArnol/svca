import numpy as np
import scipy as s
from svca.models.model1 import Model1
from svca.simulations.from_real import FromRealSimulation
from svca.models.io import *
from svca.util_functions import utils
import sys
from limix.utils.preprocess import covar_rescaling_factor_efficient
from copy import deepcopy


def run(data_dir, protein_index, output_dir, bootstrap_index,
        normalisation='standard', permute=False):
    # reading all data
    ####################################################################
    expression_file = data_dir + '/expressions.txt'
    position_file = data_dir+'/positions.txt'
    protein_names, phenotypes, X = utils.read_data(expression_file,
                                                   position_file)

    # import pdb; pdb.set_trace()
    protein_name = protein_names[protein_index, :]
    phenotype = phenotypes[:, protein_index]
    sel = range(phenotypes.shape[1])
    sel.remove(protein_index)
    kin_from = phenotypes[:, sel]
    N_samples = X.shape[0]

    # permuting cells
    if permute:
        perm = np.random.permutation(X.shape[0])
        X = X[perm, :]

    # do null simulation
    ####################################################################
    sim = FromRealSimulation(X, phenotype, kin_from)
    Y_sim = sim.simulate()

    # run model on simulated data
    ####################################################################
    # intrinsic and environmental term
    ####################################################################
    cterms = ['intrinsic', 'environmental']
    model = Model1(Y_sim, X, norm=normalisation, oos_predictions=0., cov_terms=cterms, kin_from=kin_from)
    model.reset_params()
    model.train_gp(grid_size=10)

    file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_local'
    write_variance_explained(model, output_dir, file_prefix)
    write_LL(model, output_dir, file_prefix)

    int_param = model.intrinsic_cov.getParams()
    env_param = model.environmental_cov.getParams()
    noise_param = model.noise_cov.getParams()

    ####################################################################
    # add cell-cell interactions
    ####################################################################
    model.add_cov(['interactions'])

    LL = np.Inf
    for i in range(5):
        if i == 0:
            int_bk = int_param
            env_bk = env_param
            noise_bk = noise_param
            scale_interactions = True
        else:
            int_bk = int_param * s.random.uniform(0.8, 1.2, len(int_param))
            local_bk = local_param * s.random.uniform(0.8, 1.2, len(env_param))
            noise_bk = noise_param * s.random.uniform(0.8, 1.2, len(noise_param))
            scale_interactions = False
        model.set_initCovs({'intrinsic': dir_bk,
                        'noise': noise_bk,
                        'environmental':local_bk})
        if scale_interactions:
            model.set_scale_down(['interactions'])
        else:
            model.use_scale_down = False

        model.reset_params()
        model.train_gp(grid_size=10)
        if model.gp.LML() < LL:
            LL = model.gp.LML()
            saved_params = model.gp.getParams()

    model.gp.setParams(saved_params)
    file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_interactions'
    write_variance_explained(model, output_dir, file_prefix)
    #write_r2(model, output_dir, file_prefix)
    write_LL(model, output_dir, file_prefix)
    # write_Ks(model, output_dir, file_prefix)


if __name__ == '__main__':
    data_dir = sys.argv[1]
    output_dir = sys.argv[2]
    protein_index = int(sys.argv[3])
    bootstrap_index = sys.argv[4]
    normalisation = sys.argv[5]

    #data_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/IMC_reprocessed_median/Cy1x7/'
    ## protein_index = 19
    #protein_index = 3 # 5
    #bootstrap_index = 3
    #output_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/tests/test_new_init/'
    #normalisation='quantile'

    run(data_dir, protein_index, output_dir, bootstrap_index, normalisation)
