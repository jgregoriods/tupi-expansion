from subprocess import Popen, PIPE
import multiprocessing as mp
import numpy as np


def run_model(params):
    result = Popen(["./expand", '4400', str(params[0]), str(params[1]), str(params[2]), str(params[3]), '0.0', str(params[4])],
                    stdout=PIPE).communicate()[0]
    return float(result)


param_combinations = []
for k in np.arange(0.9, 1.1, 0.1):
    for r in np.arange(0.02, 0.04, 0.01):
        for pct_migrants in np.arange(0.4, 0.5, 0.1):
            for leap_distance in range(0, 3, 2):
                for fiss_thr in np.arange(0.5, 0.7, 0.1):
                    param_combinations.append([k, r, pct_migrants, leap_distance, fiss_thr])
print(param_combinations)


for params in param_combinations:
    print(run_model(params))


"""
pool = mp.Pool(3)
res = pool.map(run_model, param_combinations)
pool.close()

print(res)
"""