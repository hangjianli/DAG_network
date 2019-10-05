import logging
import random
import numpy as np

logging.basicConfig(filename='utils.log', level=logging.DEBUG)


# This function is probably unnecessary.
def args_for_parameters(n, p, equal_var=True):
    pass


class ModelParams:
    def __init__(self,
                 n, p,
                 df_name: str,
                 b_from_data_FLAG: bool = True,
                 df_theta_name: str = 'andes',
                 df_beta_name: str = 'andes',
                 btype: str = 'continuous',
                 equal_var_FLAG: bool = True):
        self.n = n
        self.p = p
        self.b_from_data_FLAG = b_from_data_FLAG
        self.equal_var_FLAG = equal_var_FLAG

    def gen_B(self):
        pass

    def _sim_B_from_data(self):
        pass

    def gen_omg(self, seed=100):
        if self.equal_var_FLAG:
            omega = np.ones(self.p)
            omega_sq = np.ones(self.p)
            return omega, omega_sq
        np.random.seed(seed)
        # Things below are probably unnecessary
        omega = np.random.uniform(0.01, 1, self.p)
        omega[0] = 1
        omega_sq = np.square(omega)
        return omega, omega_sq



