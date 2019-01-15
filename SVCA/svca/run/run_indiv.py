import numpy as np
from svca.models.model1 import Model1
from svca.models.io import *
from svca.util_functions import utils, util_functions
import sys
from limix.utils.preprocess import covar_rescaling_factor_efficient
import argparse
import os


def run(data_dir, protein_index, output_dir,
        normalisation='quantile', permute=False):
    # reading all data
    # ------------------------------------------------------------------------
    expression_file = data_dir + '/expressions.txt'
    position_file = data_dir+'/positions.txt'
    protein_names, phenotypes, X = utils.read_data(expression_file,
                                                   position_file)

    if protein_index is None:
        for protein_index in range(len(protein_names)):
            run_indiv(protein_names, phenotypes, X, protein_index, output_dir,
                normalisation, permute)
    else:
        for p in protein_index:
            run_indiv(protein_names, phenotypes, X, p, output_dir,
                normalisation, permute)

def run_indiv(protein_names, phenotypes, X, protein_index, output_dir,
        normalisation, permute):

    protein_name = protein_names[protein_index, :]

    print('-------------------------------------------')
    print('running model for ', protein_name[0])
    print('-------------------------------------------')

    phenotype = phenotypes[:, protein_index]
    #sel = range(phenotypes.shape[1])
    sel = [i for i in range(phenotypes.shape[1]) if i != protein_index]
    #sel.remove(protein_index)
    kin_from = phenotypes[:, sel]

    N_samples = X.shape[0]

    # permuting cells
    if permute:
        perm = np.random.permutation(X.shape[0])
        X = X[perm, :]


    # intrinsic term
    # ------------------------------------------------------------------------
    cterms = ['intrinsic', 'environmental', 'interactions']
    model = Model1(phenotype, X, norm=normalisation, oos_predictions=0., cov_terms=cterms, kin_from=kin_from, cv_ix=0)
    model.reset_params()
    model.train_gp(grid_size=10)

    file_prefix = protein_name[0] + '_' + str(0) + '_interactions'
    write_variance_explained(model, output_dir, file_prefix)
    # write_LL_grid(model, output_dir, file_prefix)


if __name__ == '__main__':
    p = argparse.ArgumentParser( description='Basic run script for SVCA' )

    # I/O
    p.add_argument( '--indir',         type=str, required=True,                 help='Input data diectory containing an expression and a position file named expressions.txt and positions.txt')
    p.add_argument( '--outdir',        type=str, default=None,                  help='Output dir for text files')
    p.add_argument( '--protein_index', type=int, nargs='+', default=None,       help='Index of the protein on which to run the model, default is to run the model on all proteins')
    p.add_argument( '--normalisation', type=str, nargs='+', default="quantile", help='Type of normalisation mehtod to use for the target gene. Default is quantile normalisation. Alternative is standardisation: std')

    args = p.parse_args()

    data_dir = args.indir
    output_dir = args.outdir
    protein_index = args.protein_index
    normalisation = args.normalisation

    if output_dir is None:
        output_dir = data_dir + '/results/'
        util_functions.make_dir(output_dir)

    run(data_dir, protein_index, output_dir, normalisation)
