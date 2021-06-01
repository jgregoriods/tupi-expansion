import numpy as np
import pandas as pd

from itertools import product
from tqdm import tqdm

from model import *


SITES = pd.read_csv('sites/tupi_filtered_100.csv')


def main():
    ORIGIN = (-61.96, -10.96)
    """
    r_vals = np.arange(0.02, 0.045, 0.005)
    e_vals = np.arange(0.2, 0.45, 0.05)
    param_list = list(product(r_vals, e_vals))
    df = pd.DataFrame(columns=['r', 'eK', 'model', 'score'])
    for i in tqdm(range(len(param_list))):
        param = param_list[i]
        m1 = Model(5000, ORIGIN, param[0], param[1], False)
        m1.run()
        m1.get_score(SITES)
        df = df.append({'r': param[0], 'eK': param[1], 'model': 'null', 'score': m1.score}, ignore_index=True)

        m2 = Model(5000, ORIGIN, param[0], param[1], True)
        m2.run()
        m2.get_score(SITES)
        df = df.append({'r': param[0], 'eK': param[1], 'model': 'forest', 'score': m2.score}, ignore_index=True)
    df = df.sort_values('score')
    df.to_csv('sim_scores.csv')
    """

    m1 = Model(5000, ORIGIN, 0.025, 0.2, False)
    m1.run()
    m1.get_score(SITES)
    m1.write('res/null_model.asc')
    m1.sites['model'] = 'null'

    m2 = Model(5000, ORIGIN, 0.025, 0.2, True)
    m2.run()
    m2.get_score(SITES)
    m2.write('res/forest_model.asc')
    m2.sites['model'] = 'forest'

    df_sim_dates = pd.concat([m1.sites, m2.sites])
    df_sim_dates.to_csv('img/sim_dates.csv')


if __name__ == '__main__':
    main()
