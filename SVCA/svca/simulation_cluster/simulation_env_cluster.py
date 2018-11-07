from svca.util_functions import cluster_utils, util_functions

import glob
import os

if __name__ == '__main__':
    # --------------------------------------------------------------------------
    # To change:
    #           - analysis directory containing one directory per input image
    #           - N_prot: number of genes / proteins in your dataset
    #           - submit_cmd: the sumbission command for your cluster
    # --------------------------------------------------------------------------
    submit_cmd = 'bsub -M 800 -o tmp_log'  # To change for your cluster settings
    analysis_dir = ''
    N_prot = 26

    # --------------------------------------------------------------------------
    # Optional:
    #       -normalisation methods implemented are
    #           - quantile
    #           - std (centering and standardising)
    #       - N_sim: number of simulations >10. 
    # --------------------------------------------------------------------------
    normalisation = 'quantile'
    N_sim = 100   # > 10 the size of the env term is taken from the bootstrap index

    # list directories
    image_dirs = sorted(glob.glob(analysis_dir+'/*'))
    for image_dir in image_dirs:
        results_directory = image_dir + '/results/'
        results_directory = util_functions.make_dir(results_directory)
        for protein_ix in range(0, N_prot):
            for bootstrap_index in range(1, N_sim):
                command_line = submit_cmd + ' ' + \
                    'python ../simulations/sim_env.py ' + \
                    image_dir + ' ' + \
                    results_directory + ' ' + \
                    str(protein_ix)+ ' ' +\
                    str(bootstrap_index) + ' '+\
                    normalisation
                os.system(command_line)
