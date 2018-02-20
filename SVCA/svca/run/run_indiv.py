import numpy as np
from svca.models.model1 import Model1
from svca.models.io import *
from svca.util_functions import utils
import sys
from limix.utils.preprocess import covar_rescaling_factor_efficient


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


    # intrinsic term
    ####################################################################
    cterms = ['intrinsic']
    model = Model1(phenotype, X, norm=normalisation, oos_predictions=0.2, cov_terms=cterms, kin_from=kin_from, cv_ix=bootstrap_index)
    model.reset_params()
    model.train_gp(grid_size=3)

    file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_intrinsic'
    write_variance_explained(model, output_dir, file_prefix)
    write_pred(model, output_dir, file_prefix)
    write_LL(model, output_dir, file_prefix)
    # write_Ks(model, output_dir, file_prefix)

    # add local term
    ####################################################################
    model.add_cov(['environmental'])
    model.reset_params()
    model.train_gp(grid_size=3)

    file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_environmental'
    write_variance_explained(model, output_dir, file_prefix)
    write_pred(model, output_dir, file_prefix)
    # write_Ks(model, output_dir, file_prefix)

    write_LL_grid(model, output_dir, file_prefix)

    # add crowding term
    ####################################################################
    model.add_cov(['interactions'])
    model.reset_params()
    model.train_gp(grid_size=3)

    file_prefix = protein_name[0] + '_' + str(bootstrap_index) + '_interactions'
    write_variance_explained(model, output_dir, file_prefix)
    write_pred(model, output_dir, file_prefix)
    write_LL_grid(model, output_dir, file_prefix)
    # write_Ks(model, output_dir, file_prefix)


if __name__ == '__main__':
    # data_dir = sys.argv[1]
    # output_dir = sys.argv[3]
    # protein_index = int(sys.argv[4])
    # bootstrap_index = sys.argv[5]
    #
    # normalisation = sys.argv[6]
    #
    # try:
    #     tmp = sys.argv[7]
    #     if tmp == 'True':
    #         perm=True
    #     else:
    #         perm=False
    # except:
    #     perm=False
	# #

    data_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/IMC_paper/res/Ay10x1/'
    # protein_index = 19
    protein_index = 23  # 5
    bootstrap_index = 3
    output_dir = '/tmp/test_svca'
    normalisation='quantile'
    perm = False

    run(data_dir, protein_index, output_dir, bootstrap_index, normalisation, perm)
