from svca.util_functions import cluster_utils, util_functions

import glob
import os

if __name__ == '__main__':
    ##################################################
    # Directory of the analysis: to change, and number of genes/proteins in the
    # dataset
    ##################################################
    submit_cmd = 'bsub -M 800 -o tmp_log'  # To change for your cluster settings
    analysis_dir = '../../examples/data/IMC_example' # To change for your image directories
    N_prot = 26

    ##################################################
    # advanced parameters to change
    ##################################################
    normalisation = 'quantile'

    ##################################################
    # list directories
    image_dirs = sorted(glob.glob(analysis_dir+'/*'))
    for image_dir in image_dirs:
        results_directory = image_dir + '/results/'
        results_directory = util_functions.make_dir(results_directory)
        for protein_ix in range(0, N_prot):
            command_line = submit_cmd + ' ' + \
                'python ../run/run_indiv.py ' + \
                image_dir + ' ' + \
                results_directory + ' ' + \
                str(protein_ix)+ ' ' +\
                normalisation + ' '
            os.system(command_line)
