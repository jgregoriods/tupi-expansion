import numpy as np
import pandas as pd
import multiprocessing as mp

from itertools import product
from tqdm import tqdm

from model import *


SITES = pd.read_csv('sites/tupi_filtered_100b.csv')
ORIGIN = (-61.96, -10.96)


def test_model(params):
    m = Model(5000, ORIGIN, params[0], params[1], params[2])
    m.run()
    m.get_score(SITES)
    model_type = 'null' if not params[2] else 'forest'
    return {'r': params[0], 'eK': params[1], 'model': model_type, 'score': m.score}


def main():
    r_vals = np.arange(0.02, 0.045, 0.005)
    e_vals = np.arange(0.2, 0.45, 0.05)
    forest_vals = [False, True]
    param_list = list(product(r_vals, e_vals, forest_vals))
    df = pd.DataFrame(columns=['r', 'eK', 'model', 'score'])
    n_cores = mp.cpu_count() - 1
    pool = mp.Pool(n_cores)
    res = pool.map(test_model, param_list)
    pool.close()
    for i in res:
        df = df.append(i, ignore_index=True)
    df = df.sort_values('score')
    df.to_csv('results/sim_scores.csv')


    m1 = Model(5000, ORIGIN, 0.025, 0.3, False)
    m1.run()
    m1.get_score(SITES)
    m1.write('results/rasters/null_model.asc')
    m1.sites['model'] = 'null'

    m2 = Model(5000, ORIGIN, 0.025, 0.3, True)
    m2.run()
    m2.get_score(SITES)
    m2.write('results/rasters/forest_model.asc')
    m2.sites['model'] = 'forest'

    df_sim_dates = pd.concat([m1.sites, m2.sites])
    df_sim_dates.to_csv('results/sim_dates.csv')


if __name__ == '__main__':
    main()
