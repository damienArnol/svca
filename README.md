# Spatial Variance Components Analysis (SVCA)

## Dependencies

### Python
- numpy
- scipy
- pandas
- rpy2 (for the notebooks)
- limix (local version in this repository)

### R (for plotting)
- ggplot2
- reshape2
- gplots
- plyr
- pheatmap

### Others 
- gcc / g++

## Installation

### Installing limix

SVCA relies on a specific version of limix found in svca_limix. You should first install this package using the setup file in svca_limix.

NB: If you are already a limix user, we recommend you install svca_limix and svca in a dedicated conda environment so there is no interference between your limix versions

```
cd svca_limix
python setup.py develop
```

### Installing svca

Then install svca
```
cd ..
cd SVCA
python setup.py develop
```

## Basic usage
### Computing spatial variance signatures for single images
Running SVCA on single image and single protein can be done as illustrated in the bash script 
`SVCA/svca/run/call_run_indiv.sh`. The script calls the `run_indiv.py` script with the following inputs:
1. `data_dir='../../examples/data/IMC_example/Cy1x7/'` directory with IMC input data
2. `output_dir='./test_svca'`  the output of the analysis is saved here
3. `protein_index=23`  select the protein to be modelled
4. `normalisation='quantile'` select the normalisation method.

For the analysis of all the images and proteins we recommend to use a cluster, this is explained in the next section. 


### Computing spatial variance signatures for multiple images

*NB: For data format, look at the example in the data/IMC_example directory, which should correspond to your analysis_dir folder*

We recommend using a cluster for this.
1. Adapt the file `SVCA/svca/run_cluster/run_all_cluster.py`, to the queuing system used by your cluster.
2. Your analysis directory `analysis_dir` should contain one directory per image on which you are fitting svca
3. Each image folder should contain a `positions.txt` and an `expressions.txt`. Rows are cells and columns are (x,y) coordinates for the positions and genes for the expressions, with the gene names as the header for the expression file. No header for the positions.
4. Run `python run_all_cluster.py` in the `run_cluster` directory.
5. Results are in a `results` directory in each image directory


### Visualising variance signatures
1. Adapt the file `SVCA/svca/plot_scripts/plot_signatures.R` (bottom). `working_dir` should be your analysis directory and `plot_dir` the directory in which you want to save your plots.
2. run the file

### Cross validation

We recommend using a cluster for this. The procedure is the exact same procedure as the one for computing variance signatures, but the file used is `SVCA/svca/run_cluster/run_cross_validation_cluster.py`.

### Visualising cross validation results
1. Adapt the file `SVCA/svca/plot_scripts/plot_r2_cross_validation.R` (bottom). `working_dir` should be your analysis directory and `plot_dir` the directory in which you want to save your plots.
2. run the file `plot_r2_cross_validation.R`
