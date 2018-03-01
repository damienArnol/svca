#################################################
# To modify
#################################################
source('/home/nico/Dropbox/JRC_COMBINE/svca/SVCA/svca/plot_scripts/visu_signatures_functions.R')

working_dir = '/home/nico/Dropbox/JRC_COMBINE/svca/SVCA/examples/data/IMC_example'
#dir.create(file.path(working_dir, 'plots'))
plot_dir = paste(working_dir, 'plots', sep = '/')
plot_sig_all(working_dir, plot_dir)
