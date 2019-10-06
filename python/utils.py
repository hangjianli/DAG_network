import logging
import random
import numpy as np
import pandas as pd



logging.basicConfig(filename='utils.log', level=logging.DEBUG)

class ModelParams:
    def __init__(self,
                 n, p,
                 b_from_data_FLAG: bool = True,
                 df_theta_name: str = 'andes',
                 df_beta_name: str = 'andes',
                 btype: str = 'continuous',
                 ncopy: int = 1,
                 equal_var_FLAG: bool = True):
        self.n = n
        self.p = p
        self.b_from_data_FLAG = b_from_data_FLAG
        self.equal_var_FLAG = equal_var_FLAG
        self.df_beta_name = df_beta_name
        self.btype = btype
        self.ncopy = ncopy
        self.df_theta_name = df_theta_name

    def _construct_B_from_data(self):
        df = pd.read_csv(f'../data/BNRepo/{self.df_beta_name}.csv', index_col=0)
        pp = df.shape[0]
        if self.ncopy > 1:
            pass

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



