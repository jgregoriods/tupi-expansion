import pandas as pd

from model import *


SITES = pd.read_csv('sites/tupi_filtered.csv')


def main():
    """
    df = pd.DataFrame(columns=['r', 'eK', 'forest', 'score'])
    for r in np.arange(0.02, 0.05, 0.01):
        for e_K in np.arange(0.3, 0.5, 0.1):
            for forest in [True, False]:
                m = Model(5000, (-61.96, -10.96), r, e_K, forest)
                m.run()
                m.get_score(SITES)
                df = df.append({'r': r, 'eK': e_K, 'forest': forest, 'score': m.score},
                                ignore_index=True)
    df.to_csv('sim_scores.csv')
    """
    null_model = Model(5000, (-61.96, -10.96), 0.025, 0.3, False)
    null_model.run()
    null_model.get_score(SITES)
    print(null_model.score)
    null_model.write('res/null_model.asc')

    forest_model = Model(5000, (-61.96, -10.96), 0.025, 0.3, True)
    forest_model.run()
    forest_model.get_score(SITES)
    print(forest_model.score)
    forest_model.write('res/forest_model.asc')

if __name__ == '__main__':
    main()
