import numpy as np
from svca.util_functions import utils
import sys
from limix.utils.preprocess import covar_rescaling_factor


def write_variance_explained(trained_model, output_dir, file_prefix):
    covar_terms = trained_model.covar_terms
    var_explained = np.zeros(len(covar_terms))
    i = 0

    for cov in covar_terms.values():
        try:
            effect = 1./covar_rescaling_factor(cov.K())
        except:
            effect=0.
        var_explained[i] = effect
        i+=1

    # write file
    file_name = output_dir + '/' + file_prefix + '_effects.txt'
    result_header = ' '.join(covar_terms.keys())
    with open(file_name, 'w') as f:
        np.savetxt(f,
                   var_explained[None, :],
                   delimiter=' ',
                   header=result_header,
                   fmt='%s',
                   comments='')

def write_r2(model, output_dir, file_prefix):
    r2 = np.array(model.r2())
    r2 = np.reshape(r2, [1,1])
    r2_header = 'r2'

    file_name = output_dir + '/' + file_prefix + '_r2.txt'
    with open(file_name, 'w') as f:
        np.savetxt(f,
                   r2,
                   delimiter=' ',
                   header=r2_header,
                   fmt='%s',
                   comments='')


def write_pred(model, output_dir, file_prefix):
    pred = model.pred()
    pred_header = 'Truth Prediction'

    file_name = output_dir + '/' + file_prefix + '_pred.txt'
    with open(file_name, 'w') as f:
        np.savetxt(f,
                   pred,
                   delimiter=' ',
                   header=pred_header,
                   fmt='%s',
                   comments='')

def write_LL(model, output_dir, file_prefix):
    try:
        LL = np.array(-model.gp.LML())
    except:
        LL=np.array([np.nan])
    LL = np.reshape(LL, [1,1])
    LL_header = 'LL'

    file_name = output_dir + '/' + file_prefix + '_LL.txt'
    with open(file_name, 'w') as f:
        np.savetxt(f,
                   LL,
                   delimiter=' ',
                   header=LL_header,
                   fmt='%s',
                   comments='')

def write_LL_grid(model, output_dir, file_prefix):
    LML = model.LMLs
    g = model.l_grid

    res = np.concatenate((g[:,None], LML[:,None]), axis=1)
    LL_header = 'length LML'

    file_name = output_dir + '/' + file_prefix + '_LL_grid.txt'
    with open(file_name, 'w') as f:
        np.savetxt(f,
                   res,
                   delimiter=' ',
                   header=LL_header,
                   fmt='%s',
                   comments='')
