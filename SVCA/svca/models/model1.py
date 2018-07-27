from model_base import Model

from svca.util_functions import utils
import utils_loc
import numpy as np

# import useful limix functions
from limix.core.covar import SQExpCov
from limix.core.covar import ZKZCov
from limix.core.covar import FixedCov
from limix.core.covar import SumCov

import limix.core.covar.combinators

from limix.core.mean import MeanBase

from limix.core.gp import GP

from limix.utils.preprocess import covar_rescaling_factor_efficient


class Model1(Model):
    """docstring for Model1."""
    def __init__(self, Y, positions, types=None, norm='quantile', oos_predictions=0., cov_terms=None, kin_from=None, cv_ix=None):
        super(Model1, self).__init__()

        self.init_covs= dict()
        self.init_from_previous = False
        self.use_scale_down = False
        self.scale_down = None

        self.Y = Y
        self.n_samples = len(Y)
        self.types = types
        self.X = positions
        self.norm = norm
        self.kin_from=kin_from

        if cv_ix is not None:
            self.cv_ix = float(int(cv_ix))
        else:
            self.cv_ix = cv_ix

        self.def_training_set(oos_predictions)

        if cov_terms is None:
            cov_terms = ['intrinsic', 'interactions', 'environmental']

        self.init_model(cov_terms)


    ##########################
    # Building Model
    ##########################
    def build_Kinship(self):
        if self.kin_from is None:
            assert self.types is not None, 'Provide a vector of cell types or a cell state matrix to compute kinship'
            self.Kin = utils.build_cell_type_kinship(self.types)
        else:
            self.kin_from -= self.kin_from.mean(axis=0)
            #self.kin_from /= self.kin_from.std(axis=0)
            self.Kin = self.kin_from.dot(self.kin_from.transpose())
            self.Kin *= covar_rescaling_factor_efficient(self.Kin)

    def build_cov(self, cov_terms):
        # import pdb; pdb.set_trace()
        self.intrinsic_cov = self.build_intrinsic()
        self.interactions_cov = self.build_interactions()
        self.environmental_cov = self.build_environmental()
        self.noise_cov = self.build_noise()

        self.covar_terms = dict()
        if 'intrinsic' in cov_terms:
            self.covar_terms['intrinsic'] = self.intrinsic_cov
        if 'interactions' in cov_terms:
            self.covar_terms['interactions'] = self.interactions_cov
        if 'environmental' in cov_terms:
            self.covar_terms['environmental'] = self.environmental_cov

        self.covar_terms['noise'] = self.noise_cov
        self.covar = apply(SumCov, self.covar_terms.values())
        self.reset_params()

    def build_intrinsic(self):
        K = self.Kin[self.train_set,:][:,self.train_set]
        K_cross = self.Kin[~self.train_set, :][:, self.train_set]
        if K_cross.shape == (0,0):
            K_cross = None
        return FixedCov(K, K_cross)

    def build_interactions(self):
        # build zkz kernel
        Xstar = ~self.train_set
        if not any(Xstar): Xstar = None

        ###############
        # toto = np.ones([self.n_samples, self.n_samples])
        # interactions_tmp = ZKZCov(self.X, toto, Xstar=Xstar)
        interactions_tmp = ZKZCov(self.X, self.Kin, Xstar=Xstar)
        #####################

        interactions_tmp.act_scale = True
        interactions_tmp.act_length = False
        interactions_tmp.length =15.

        return interactions_tmp

    def build_environmental(self):
        Xtrain, Xstar = self.X[self.train_set, :], self.X[~self.train_set, :]
        if Xstar.shape == (0,0): Xstar = None
        environmental_cov = SQExpCov(Xtrain, Xstar=Xstar)
        environmental_cov.act_length = False

        return environmental_cov

    def build_noise(self):
        # TODO one per type ?
        n_train = self.train_set.sum()
        glob_cov = FixedCov(np.eye(n_train))
        return glob_cov

    def add_cov(self, cov_term):
        # import pdb; pdb.set_trace()
        for term in cov_term:
            if term in self.covar_terms:
                continue
            elif term == 'intrinsic':
                self.covar_terms['intrinsic'] = self.intrinsic_cov
            elif term == 'interactions':
                self.covar_terms['interactions'] = self.interactions_cov
            elif term == 'environmental':
                self.covar_terms['environmental'] = self.environmental_cov
            else:
                raise Exception('covariance term not recognised ')

        self.covar = apply(SumCov, self.covar_terms.values())
        self.build_gp()

    def rm_cov(self, cov_term):
        for term in cov_term:
            if term in self.covar_terms:
                del self.covar_terms[term]
            else:
                print 'cov term ', term, ' not found for deletion '
        self.covar = apply(SumCov, self.covar_terms.values())
        self.build_gp()

    def set_initCovs(self, cov_dir):
        # input a dictionary
        self.init_covs = cov_dir
        self.init_from_previous = True


    def reset_from_previous(self):
        for cov_key in self.covar_terms.keys():
            if cov_key in self.init_covs.keys():
                self.covar_terms[cov_key].setParams(self.init_covs[cov_key])

            else:
                self.covar_terms[cov_key].scale = 1.
                k = covar_rescaling_factor_efficient(self.covar_terms[cov_key].K())
                self.covar_terms[cov_key].scale = 0.5 * k

        if self.use_scale_down:
            for term in self.scale_down:
                self.covar_terms[term].scale *= 1e-10

        self.covar = apply(SumCov, self.covar_terms.values())
        self.build_gp()

    def set_scale_down(self, terms):
        self.scale_down = terms
        self.use_scale_down = True

    def reset_params(self):

        if self.init_from_previous:
            self.reset_from_previous()
            return

        # import pdb; pdb.set_trace()
        self.intrinsic_cov.scale = 1.

        self.interactions_cov.scale = 1.

        self.environmental_cov.scale = 1.
        self.noise_cov.scale = 1.

        used_covar = self.covar_terms.values()
        n = len(used_covar)

        for cov in used_covar:
            k = covar_rescaling_factor_efficient(cov.K())/n
            # TODO: not so clean, would be good to have setInterParams ?
            new_params = np.log(cov.getInterParams() * k)
            cov.setParams(new_params)

        if self.use_scale_down:
            for term in self.scale_down:
                self.covar_terms[term].scale *= 1e-6

    ##########################
    # Train model
    ##########################
    # def train_gp(self):
        # self.gp.optimize()

    def train_gp(self, grid_size=20):
        # if no interactions  in covar ters, no need for grid search
        g = ('interactions' in self.covar_terms or 'environmental' in self.covar_terms)
        if not g:
            self.l_grid = None
            self.LMLs = None
            try:
                self.gp.optimize()
            except:
                # import pdb; pdb.set_trace()
                print 'No convergence for non-grid optimisation in Model1'
                exit()
            return

        #itterate length scale of interactions term (shared parameter)
        self.l_grid = utils_loc.get_l_grid(self.X, grid_size)
        self.LMLs = np.zeros(grid_size)

        # NOTE reducing grid for debugging
        # self.l_grid = self.l_grid[[2]]
        # self.l_grid = np.array([20.])
        # self.LMLs = np.zeros([1])

        LML = np.Inf
        params_saved = None
        best_l = 0.
        i = 0

        for l in self.l_grid:
            self.interactions_cov.length = l**2.
            self.environmental_cov.length = l**2.

            self.reset_params()

            try:
                self.gp.optimize()
            except:
                # import pdb; pdb.set_trace()
                self.LMLs[i] = np.nan
                i+=1
                continue
            self.LMLs[i] = self.gp.LML()
            i+=1
            if self.gp.LML() < LML:
                best_l = l
                LML = self.gp.LML()
                params_saved = self.gp.getParams()

        self.interactions_cov.length = best_l**2.
        self.environmental_cov.length = best_l**2.
        try:
            self.gp.setParams(params_saved)
        except:
            exit(1)
