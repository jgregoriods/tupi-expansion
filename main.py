import numpy as np
import pandas as pd
import multiprocessing as mp

from itertools import product

from model import *


SITES = pd.read_csv('sites/tupi_dates.csv')
ORIGIN = (-61.96, -10.96)


def test_model(params):
    """
    Runs simulations with different parameter combinations and calculates their
    score.
    """
    model = Model(5000, ORIGIN, params[0], params[1], params[2])
    model.run()
    model.get_score(SITES)
    model_type = 'null' if not params[2] else 'forest'
    return {'growth': params[0], 'emigration': params[1], 'model': model_type, 'score': model.score}


def main():
    growth_vals = np.arange(0.02, 0.045, 0.005)
    emigration_vals = np.arange(0.2, 0.45, 0.05)
    forest_vals = [False, True]
    param_list = list(product(growth_vals, emigration_vals, forest_vals))
    df = pd.DataFrame(columns=['growth', 'emigration', 'model', 'score'])
    n_cores = mp.cpu_count() - 1
    pool = mp.Pool(n_cores)
    print(f'Running {len(param_list)} simulations. This may take a while ...')
    scores = pool.map(test_model, param_list)
    pool.close()
    for i in scores:
        df = df.append(i, ignore_index=True)
    df = df.sort_values('score')
    df.to_csv('results/sim_scores.csv')

    null_model = Model(5000, ORIGIN, 0.025, 0.3, False)
    null_model.run()
    null_model.get_score(SITES)
    null_model.write('results/rasters/null_model.asc')
    null_model.sites['model'] = 'null'

    forest_model = Model(5000, ORIGIN, 0.025, 0.3, True)
    forest_model.run()
    forest_model.get_score(SITES)
    forest_model.write('results/rasters/forest_model.asc')
    forest_model.sites['model'] = 'forest'

    df_sim_dates = pd.concat([null_model.sites, forest_model.sites])
    df_sim_dates.to_csv('results/sim_dates.csv')


if __name__ == '__main__':
    main()
