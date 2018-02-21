from svca.util_functions import cluster_utils, util_functions

import glob
import os

if __name__ == '__main__':

    ##################################################
    # Directory of the analysis: to change, and number of genes/proteins in the
    # dataset
    ##################################################
    analysis_dir = '/gpfs/nobackup/stegle/users/arnol/spatial/simulations/IMC_env_simulations/'
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
            bootstrap_index = 1
            command_line = \
                'bsub -q research-rh7 -o tmp_log -M 800 -R "rusage[mem=800]" python ../run/run_indiv.py ' + \
                image_dir + ' ' + \
                results_directory + ' ' + \
                str(protein_ix)+ ' ' +\
                str(bootstrap_index) + ' '+\
                normalisation + ' '
            os.system(command_line)
