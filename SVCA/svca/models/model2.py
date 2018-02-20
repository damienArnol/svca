from model1 import Model1

from svca.util_functions import utils
import utils_loc
import numpy as np

# import useful limix functions
from limix.core.covar import SQExpCov
from limix.core.covar import ZKZCov
from limix.core.covar import FixedCov
from limix.core.covar import SumCov
from limix.core.covar import ProdCov
from limix.core.covar import CategoricalCov
from limix.core.covar import ProdCov

from limix.core.mean import MeanBase

from limix.core.gp import GP

from limix.utils.preprocess import covar_rescaling_factor_efficient

from copy import deepcopy



class Model2(Model1):
    """docstring for Model1."""
    def __init__(self, Y, types, positions, norm='quantile', oos_predictions=0., cov_terms=None, kin_from=None, cv_ix=None):


        super(Model2, self).__init__(Y, types, positions, norm, oos_predictions, cov_terms, None, kin_from, cv_ix)









        # import pdb; pdb.set_trace()
        # self.loc_noise_cov.length = 15.
        # self.loc_noise_cov.setPenalty(15., 1.)


    ##########################
    # Train model
    ##########################
    # def train_gp(self):
        # self.gp.optimize()
