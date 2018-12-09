from svca.util_functions import utils
import svca.models.utils_loc as utils_loc
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

from limix.utils.preprocess import covar_rescaling_factor


class Model(object):
    """
        Model is a general class for building and training a spatial variance model.
        Contains all the functions which are not specific to a given model
    """
    def __init__(self):
        pass

    ##########################
    # Preprocessing steps
    ##########################
    '''
    Normalisation of Y
    '''
    def preprocess_input(self):
        # normalise phenotype
        if self.norm == 'quantile':
            # import pdb; pdb.set_trace()
            self.Y = utils.quantile_normalise_phenotype(self.Y)
        elif self.norm == 'std':
            self.Y = utils.normalise_phenotype(self.Y)
        else:
            raise Exception('normalisation method not understood')

    '''
    Define a training set and a test set for out of sample prediction
    '''
    def def_training_set(self, oos_predictions):
        if self.cv_ix is None:
            tmp = np.array([True for i in range(self.n_samples)])
            if oos_predictions == 0.:
                self.train_set = tmp
            elif 0. < oos_predictions < 1.:
                test_ix = np.random.choice(range(self.n_samples), int(oos_predictions * self.n_samples), replace=False)
                tmp[test_ix] = False
                self.train_set = tmp
            else:
                raise Exception('oos_predictions out of range, should be in [0;1[')

        else:
            # set seed and get an index permutation and step size
            np.random.seed(0)
            permuted_indices = np.random.permutation(self.X.shape[0])
            step_size = len(permuted_indices) * oos_predictions

            # select test set
            first_ix = int(self.cv_ix * step_size)
            last_ix = int(self.cv_ix  * step_size + step_size)
            test_set = permuted_indices[first_ix:last_ix]

            # define boolean vector for train set
            self.train_set = np.array([True for i in range(self.n_samples)])
            self.train_set[test_set] = False

    ##########################
    # Building Model
    ##########################
    '''
        General way of initialising a mdel:
    '''
    def init_model(self, cov_terms):
        self.preprocess_input()   # defined in parent
        self.build_Kinship()
        self.build_cov(cov_terms)
        self.build_mean()
        self.build_gp()

    '''
        The following functions are specific to a given model and have to be implemented
        in the relevant children class
    '''
    def build_Kinship(self):
        pass

    def build_cov(self):
        pass

    def add_cov(self):
        pass

    def rm_cov(self):
        pass

    '''
        General way to build the mean term of a GP model for limix
    '''
    def build_mean(self):
        Y_tmp = self.Y
        Y_tmp = Y_tmp[self.train_set, :]
        self.mean = MeanBase(Y_tmp)

    '''
        Creating a limix GP object
    '''
    def build_gp(self):
        self.gp = GP(self.mean, self.covar)

    ##########################
    # Train model
    ##########################
    '''
        The way the model is trained is specific to the model and has to be implemented
        in the relevant classes
    '''
    def train_gp(self):
        pass

    ##########################
    # Prediction from model
    ##########################
    '''
        General functions for out of sample prediction
    '''
    def predict(self):
        try:
            return self.gp.predict()
        except:
            return np.array([[np.nan]])

    def r2(self):
        Y_pred = self.predict()[:,0]
        Y_true = self.Y[:,0][~self.train_set]

        res = ((Y_true - Y_pred)**2.).sum()
        var = ((Y_true - Y_true.mean())**2.).sum()

        return 1. - res/var

    def pred(self):
        Y_pred = self.predict()[:,0]
        Y_true = self.Y[:,0][~self.train_set]

        return np.concatenate((Y_true[:, None], Y_pred[:, None]), axis=1)
