from svca.util_functions import cluster_utils, util_functions

import glob
import os

if __name__ == '__main__':
    ##################################################
    # Directory of the analysis: to change, and number of genes/proteins in the
    # dataset
    ##################################################
    analysis_dir = '/gpfs/nobackup/stegle/users/arnol/spatial/tests_rep_paper/IMC_data/'
    N_prot = 26

    ##################################################
    # advanced parameters to change
    ##################################################
    normalisation = 'quantile'
    N_sim = 100

    # list directories
    image_dirs = sorted(glob.glob(analysis_dir+'/*'))
    for image_dir in image_dirs:
        results_directory = image_dir + '/results/'
        results_directory = util_functions.make_dir(results_directory)
        for protein_ix in range(0, N_prot):
            for bootstrap_index in range(1, N_sim):
                command_line = \
                    'bsub -q research-rh7 -o tmp_log -M 800 -R "rusage[mem=800]" python ../simulations/sim_null.py ' + \
                    image_dir + ' ' + \
                    results_directory + ' ' + \
                    str(protein_ix)+ ' ' +\
                    str(bootstrap_index) + ' '+\
                    str(N_sim) + ' '+\
                    normalisation + ' '
                os.system(command_line)
