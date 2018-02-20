import numpy as np

from svca.models.model2 import Model2
# import useful limix functions
from limix.utils.preprocess import covar_rescaling_factor_efficient
from limix.core.covar import ZKZCov


class FromRealSimulation(object):
    def __init__(self, X, Y, kin_from):
        self.X = X
        self.Y = Y
        self.kin_from = kin_from

    def simulate(self, cterms=['direct', 'local'], env_size=None): # Should we add local ?
        # train gp with requested terms
        model = Model2(self.Y, ['all'], self.X, norm='quantile', oos_predictions=0.,
                       cov_terms=cterms, kin_from=self.kin_from)
        model.reset_params()
        model.train_gp(grid_size=10)

        # simulate from gp
        k = model.gp.covar.K()
        k *= covar_rescaling_factor_efficient(k)

        # manually add a cross-talk term
        if env_size is not None:
            assert 0. < env_size < 1., 'env size must be between 0 and 1 '
            zkz = ZKZCov(self.X, model.Kin)
            zkz.length = model.loc_noise_cov.length
            tmp = zkz.K()
            tmp *= covar_rescaling_factor_efficient(tmp)

            tmp *= (env_size / (1- env_size))
            k += tmp

        res = np.random.multivariate_normal([0.]*k.shape[0], k)

        return res

class FromRealSimulation3(object):
    def __init__(self, X, Y, kin_from):
        self.X = X
        self.Y = Y
        self.kin_from = kin_from

    def simulate(self, cterms=['direct', 'local', 'env'], env_size=None):
        # train gp with requested terms
        model = Model2(self.Y, ['all'], self.X, norm='quantile', oos_predictions=0.,
                       cov_terms=cterms, kin_from=self.kin_from)
        model.reset_params()
        model.train_gp(grid_size=10)

        # simulate from gp after removing env term
        k = model.covar_terms['direct'].K() + \
            model.covar_terms['local'].K() + \
            model.covar_terms['noise'].K()
        k *= covar_rescaling_factor_efficient(k)

        # manually add a cross-talk term
        if env_size is not None:
            assert 0. < env_size < 1., 'env size must be between 0 and 1 '
            tmp = model.covar_terms['env'].K()
            tmp *= covar_rescaling_factor_efficient(tmp)

            tmp *= (env_size / (1. - env_size))
            k += tmp

        res = np.random.multivariate_normal([0.]*k.shape[0], k)

        return res


class FromRealSimulation2(object):
    def __init__(self, X, Y, kin_from):
        self.X = X
        self.Y = Y

        self.kin_from = kin_from

    def simulate(self, cterms=['direct', 'local'], env_size=None): # Should we add local ?
        # train gp with requested terms
        model = Model2(self.Y, ['all'], self.X, norm='quantile', oos_predictions=0.,
                       cov_terms=cterms, kin_from=self.kin_from)

        # do not train model but take variance explained evenly distribuuted
        # ac
        model.env_cov.length =100
        model.crowd_cov.length = 100
        model.loc_noise_cov.length = 100

        model.reset_params()

        # simulate from gp
        k = model.gp.covar.K()
        k *= covar_rescaling_factor_efficient(k)

        # manually add a cross-talk term
        if env_size is not None:
            assert 0. < env_size < 1., 'env size must be between 0 and 1 '
            zkz = ZKZCov(self.X, model.Kin)
            zkz.length = model.loc_noise_cov.length
            tmp = zkz.K()
            tmp *= covar_rescaling_factor_efficient(tmp)

            tmp *= (env_size / (1- env_size))
            k += tmp

        res = np.random.multivariate_normal([0.]*k.shape[0], k)

        return res


if __name__ == '__main__':
    x = np.reshape(np.random.randn(20), [10,2])
    y = np.random.randn(10)
    kin_from = np.reshape(np.random.randn(70), [10,7])

    tmp = FromRealSimulation(x, y, kin_from)
    toto = tmp.simulate()
    print toto
