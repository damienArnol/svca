
# from spatialGP.util_functions import utils
# import utils_loc
import numpy as np
import scipy.spatial as ss

# import useful limix functions
from limix.core.covar import SQExpCov
from limix.core.covar import ZKZCov
from limix.utils.preprocess import covar_rescaling_factor_efficient

class Simulation(object):
    """docstring for Simulation."""
    def __init__(self,
                 n_samples,
                 n_genes=10,
                 l2=100,
                 exp=None,
                 X=None):

        self.n_samples = n_samples
        self.n_genes = n_genes
        self.l2 = l2

        # intialise predictors
        self.initialise(exp, X)

    def initialise(self, exp, X):
        # initialise positions
        if X is None:
            L = np.sqrt(self.n_samples/0.004)
            X1 = np.random.uniform(0, L, [self.n_samples/2,2])  # problem is crowding is uniform
            X2 = np.reshape(L * np.random.beta(1, 7, self.n_samples), [self.n_samples/2, 2])
            X = np.concatenate((X1, X2), axis =0)
        self.X = X
        self.n_samples = self.X.shape[0]


        # compute squared distance between cells and set diagonal to infinity
        rv = ss.distance.pdist(self.X,'euclidean')**2.
        self.d2 = ss.distance.squareform(rv)
        self.d2[self.d2==0] = np.Inf

        # define expression of input genes
        if exp is None:
            exp = np.random.randn(self.n_samples*self.n_genes)
            exp = np.reshape(exp, [self.n_samples,self.n_genes])
        exp -= exp.mean(axis=0)
        # exp /= exp.std(axis=0)
        self.exp = exp
        self.n_genes = exp.shape[1]


class MechaSimulation(Simulation):
    def __init__(self,
                 n_samples,
                 n_genes=10,
                 l2=100,
                 exp=None,
                 X=None):

        super(MechaSimulation, self).__init__(n_samples, n_genes, l2, exp, X)

    # NOTE for now only env and spill size matter
    def simulate(self,
                 env_size=0.5,
                 loc_size=0.,
                 int_size=1.,
                 spillover_size=0.5,
                 noise_size=1.,
                 c_c_number=5.):

        assert env_size <= 1., 'choose env size smaller than 1'
        assert spillover_size <= 1., 'choose spillover size smaller than 1'

        self.env_size = env_size
        self.loc_size = loc_size
        self.intrinsic_size = int_size
        self.spillover_size = spillover_size
        self.noise_size = noise_size
        self.c_c_number = c_c_number

        self.Y = np.zeros(self.n_samples)

        self.simulate_intrinsic()
        # self.simulate_local()
        self.simulate_env()

        # add noise
        self.noise_effect = np.random.randn(self.n_samples)
        self.noise_effect /= self.noise_effect.std()

        # self.Y_truth = self.Y.copy()
        # TODO check that this is ok for %tage of variance explained (eg check r2)
        res = self.int_effect + self.noise_effect
        res_std = res.std()
        self.int_effect /= res_std
        self.noise_effect /= res_std
        self.Y = np.sqrt(env_size) * self.env_effect + np.sqrt(1 - self.env_size) * (self.int_effect + self.noise_effect)

        # ad spilover effect
        self.simulate_spillover()
        self.Y = np.sqrt(1 - self.spillover_size) * self.Y + np.sqrt(self.spillover_size) *  self.spill_effect

        #
        # # readd noise
        # self.noise_effect2 = np.random.randn(self.n_samples)
        # self.Y += self.noise_effect2

    #
    def simulate_spillover(self):
        assert self.spillover_size <= 1., 'choose spillover size between 0 and 1'

        # build probability matrix for spillover
        prob = 1. / self.d2
        prob /= prob.sum(axis=1)[:,None]
        # np.random.choice(range(len(prob)), p=prob, replace=False, size=2)

        # select cells for spillover
        cell_sel = [np.random.choice(range(prob.shape[0]), p=prob[i,:], replace=False, size=2) for i in range(prob.shape[0])]
        cell_sel_bin = np.zeros([self.exp.shape[0], self.exp.shape[0]])
        for i in range(self.exp.shape[0]):
            cell_sel_bin[i, cell_sel[i]] = 1.
        cell_sel_bin = cell_sel_bin + cell_sel_bin.T
        cell_sel_bin[cell_sel_bin != 0] = 1.
        cell_sel_bin /= cell_sel_bin.sum(axis=1)[:,None]

        # averaging effect between selected cells for X
        self.spill_effect_exp = cell_sel_bin.dot(self.exp)
        self.exp = np.sqrt(1 - self.spillover_size) * self.exp + np.sqrt(self.spillover_size) * self.spill_effect_exp

        # averaging effect between selected cells for Y
        self.spill_effect = cell_sel_bin.dot(self.Y)
        self.spill_effect /= self.spill_effect.std()

    def simulate_intrinsic(self):
        # kin = self.exp.dot(self.exp.transpose())
        # self.Y += self.intrinsic_size * np.random.multivariate_normal([0.]*self.n_samples, kin)
        self.int_effect = self.exp.dot(np.random.randn(self.n_genes))
        self.int_effect /= self.int_effect.std()

    def simulate_env(self):
        # use a function different from SE
        tmp = self.d2.copy()

        tmp=self.select_min_d(tmp, self.c_c_number)
        f = 1./tmp

        # use SE
        # se = SQExpCov(self.X)
        # se.length = self.l2
        # k = se.K()
        # k -= k.diagonal() * np.eye(k.shape[0])
        # eff = k.dot(self.exp)
        eff = f.dot(self.exp)
        eff -= eff.mean(axis =0)
        eff /= eff.std(axis =0)
        # tmp = f.dot(self.exp.dot(self.effects))

        self.env_effect = eff.dot(np.random.randn(self.n_genes))
        self.env_effect /= self.env_effect.std()

    def select_min_d(self, d2, n):
        distance_threshold = np.sort(d2, axis = 1)[:, n-1]
        sim_dis_f = d2.copy()
        for i in range(d2.shape[0]):
            sim_dis_f[i,sim_dis_f[i,:]> distance_threshold[i]] = np.inf
            sim_dis_f[i,i] = np.inf
        return(sim_dis_f)


class GPSimulation(Simulation):
    def __init__(self,
                 n_samples,
                 n_genes=10,
                 env_size=2.,
                 l2=100):
        super(GPSimulation, self).__init__(n_samples, n_genes, env_size, l2)

    def simulate(self, terms='all'):
        if terms == 'all':
            self.simul_terms = ['intrinsic', 'crowding', 'local', 'env']
        else:
            self.simul_terms = terms

        self.covar = np.zeros([self.n_samples, self.n_samples])

        if 'intrinsic' in self.simul_terms:
            self.simulate_intrinsic()
        if 'crowding' in self.simul_terms:
            self.simulate_crowding()
        if 'local' in self.simul_terms:
            self.simulate_local()
        if 'env' in self.simul_terms:
            self.simulate_env()

        # add noise and simulate
        self.covar += np.eye(self.n_samples)
        self.Y = np.random.multivariate_normal([0.]*self.n_samples, self.covar)

    def simulate_intrinsic(self):
        kin = self.exp.dot(self.exp.transpose())
        kin *= covar_rescaling_factor_efficient(kin)
        self.covar += kin

    def simulate_crowding(self):
        k1 = np.ones([self.n_samples, self.n_samples])
        # kin *= covar_rescaling_factor_efficient(kin)

        tmp = ZKZCov(self.X, k1)
        tmp.length = self.l2
        k = tmp.K()
        k *= covar_rescaling_factor_efficient(k)

        self.covar += k

    def simulate_local(self):
        tmp = SQExpCov(self.X)
        tmp.length = self.l2
        k = tmp.K()
        k *= covar_rescaling_factor_efficient(k)

        self.covar += k

    def simulate_env(self):
        kin = self.exp.dot(self.exp.transpose())
        kin *= covar_rescaling_factor_efficient(kin)

        tmp = ZKZCov(self.X, kin)
        tmp.length = self.l2
        k = tmp.K()
        k *= covar_rescaling_factor_efficient(k)

        self.covar += k



if __name__ == '__main__':
    data_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/IMC_reprocessed_median/Cy1x7/'
    # # protein_index = 19
    # protein_index = 3 # 5
    # cell_types_file = ''
    # env_size = 3
    # output_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/tests/test_new_init/'
    # normalisation='quantile'
    # perm = False
    from spatialGP.simulations.simulation_models import *
    exp = np.loadtxt(data_dir+'/expressions.txt', skiprows=1)
    X = np.loadtxt(data_dir+'/positions.txt', delimiter=',')

    sim = MechaSimulation(n_samples=0, n_genes=0, env_size=0., l2=100, exp=exp, X=X)
    sim.simulate()

    # run svca
    from spatialGP.models.model2 import Model2
    cterms = ['direct', 'local', 'env']
    model = Model2(sim.Y, None, X, norm='quantile', oos_predictions=0., cov_terms=cterms, kin_from=X)
