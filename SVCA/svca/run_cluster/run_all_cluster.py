from svca.util_functions import cluster_utils, util_functions

import glob
import os

if __name__ == '__main__':
    # all_expressions_file = '/homes/arnol/random_effect_model/data/membrane_expressions'
    # IMC_dir='/hps/nobackup/stegle/users/arnol/data/SampleSet_1/IMC_data/'
    #expression_dir = '/homes/arnol/current/data_selection/data/cytoplasm_expression/'
    #positions_dir = '/homes/arnol/current/data_selection/data/positions/'
    analysis_dir = '/gpfs/nobackup/stegle/users/arnol/spatial/simulations/IMC_env_simulations/'
    #cell_types_dir = '/homes/arnol/current/data_selection/data/cell_types/'

    #split_expression_file(all_expressions_file, IMC_dir, expression_dir)
    #cluster_utils.create_analysis_tree(positions_dir, expression_dir, analysis_dir)

    # parameters to change
    normalisation = 'quantile'
    N_bootstraps = 100
    rank=1
    permutation=False
    N_prot = 26  # 140

    # list directories
    image_dirs = sorted(glob.glob(analysis_dir+'/*'))
    image_dirs = image_dirs[::5]
    #type_files = sorted(glob.glob(cell_types_dir+'/*'))
    for image_dir in image_dirs:
        type_file = 'blank'
        results_directory = image_dir + '/results/'
        results_directory = util_functions.make_dir(results_directory)
        for protein_ix in range(0, N_prot):
            for bootstrap_index in range(1, N_bootstraps):
                command_line = \
                    'bsub -q research-rh7 -o tmp_log -M 800 -R "rusage[mem=800]" python ../simulations/sim_from_real.py ' + \
                    image_dir + ' ' + \
                    type_file + ' ' +\
                    results_directory + ' ' + \
                    str(protein_ix)+ ' ' +\
                    str(bootstrap_index) + ' '+\
                    normalisation + ' '+\
                    str(permutation)
                os.system(command_line)
