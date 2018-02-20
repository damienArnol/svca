import numpy as np
import scipy as s
from svca.models.model2 import Model2
from svca.simulations.from_real import FromRealSimulation3
from svca.models.io import *
from svca.util_functions import utils
import sys
from limix.utils.preprocess import covar_rescaling_factor_efficient
from copy import deepcopy


def run(data_dir, protein_index, cell_types_file, output_dir, bootstrap_index,
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

    # cell_types, cell_types_names, n_types, types_by_name = utils.read_types(cell_types_file, with_names=True)

    # X, phenotype, cell_types, types_by_name =  utils.subset_image(X, phenotype, cell_types, types_by_name, order=True, size=10000)
    N_samples = X.shape[0]

    # permuting cells
    if permute:
        perm = np.random.permutation(X.shape[0])
        X = X[perm, :]

    # do null simulation
    ####################################################################
    sim = FromRealSimulation3(X, phenotype, kin_from)
    Y_sim = sim.simulate()

    # run model on simulated data
    ####################################################################
    # intrinsic term
    ####################################################################
    cterms = ['direct', 'local']
    cell_types = ['all']
    model = Model2(Y_sim, cell_types, X, norm=normalisation, oos_predictions=0., cov_terms=cterms, kin_from=kin_from)
    model.reset_params()
    model.train_gp(grid_size=10)

    file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_local'
    write_variance_explained(model, output_dir, file_prefix)
    write_LL(model, output_dir, file_prefix)

    # add local term
    ####################################################################
    # model.add_cov(['local'])
    # model.reset_params()
    # model.train_gp(grid_size=10)
    #
    # file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_local'
    # write_variance_explained(model, output_dir, file_prefix)
    # write_r2(model, output_dir, file_prefix)
    # # write_Ks(model, output_dir, file_prefix)
    #
    # write_LL_grid(model, output_dir, file_prefix)

    # add crowding term
    ####################################################################
    dir_param = model.dir_cov.getParams()
    local_param = model.loc_noise_cov.getParams()
    noise_param = model.glob_noise_cov.getParams()

    model.add_cov(['env'])

    LL = np.Inf
    for i in range(5):
        if i == 0:
            dir_bk = dir_param
            local_bk = local_param
            noise_bk = noise_param
            scale_env = True
        else:
            dir_bk = dir_param * s.random.uniform(0.8, 1.2, len(dir_param))
            local_bk = local_param * s.random.uniform(0.8, 1.2, len(local_param))
            noise_bk = noise_param * s.random.uniform(0.8, 1.2, len(noise_param))
            scale_env = False
        model.set_initCovs({'direct': dir_bk,
                        'noise': noise_bk,
                        'local':local_bk})
        if scale_env:
            model.set_scale_down(['env'])
        else:
            model.use_scale_down = False

        model.reset_params()
        model.train_gp(grid_size=10)
        if model.gp.LML() < LL:
            LL = model.gp.LML()
            saved_params = model.gp.getParams()

    model.gp.setParams(saved_params)
    file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_env'
    write_variance_explained(model, output_dir, file_prefix)
    #write_r2(model, output_dir, file_prefix)
    write_LL(model, output_dir, file_prefix)
    # write_Ks(model, output_dir, file_prefix)


if __name__ == '__main__':
    data_dir = sys.argv[1]
    cell_types_file = sys.argv[2]
    output_dir = sys.argv[3]
    protein_index = int(sys.argv[4])
    bootstrap_index = sys.argv[5]

    normalisation = sys.argv[6]

    try:
        tmp = sys.argv[7]
        if tmp == 'True':
            perm=True
        else:
            perm=False
    except:
        perm=False
      #

    #data_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/IMC_reprocessed_median/Cy1x7/'
    ## protein_index = 19
    #protein_index = 3 # 5
    #cell_types_file = ''
    #bootstrap_index = 3
    #output_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/tests/test_new_init/'
    #normalisation='quantile'
    #perm = False

    run(data_dir, protein_index, cell_types_file, output_dir, bootstrap_index, normalisation, perm)
