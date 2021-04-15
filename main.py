import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set()

from mpl_toolkits.axes_grid1 import make_axes_locatable
from model import *


SITES = pd.read_csv('sites/tupi.csv')
ele = np.loadtxt('layers/ele.asc', skiprows=6)
ele[ele < 1] = np.nan
ele[ele >= 1] = 1


def plot_graphs(scores):
    fig, axes = plt.subplots(1, 2)
    for score, ax, c in zip(scores, axes.flat, ['mediumblue', 'darkred']):
        ax.scatter(SITES['dist'], SITES['bp'], c='darkgray')
        ax.scatter(score['dist'][score['sim_dates'] != 0],
                   score['sim_dates'][score['sim_dates'] != 0],
                   c=c)
        ax.set_xlabel('distance (km)')
        ax.set_ylabel('age (cal BP)')
    fig.set_size_inches(12, 6)
    fig.tight_layout()
    plt.savefig('img/plots.jpeg', dpi=300)


def plot_time_slices(slices):
    num_slices = len(slices) // 2
    fig, axes = plt.subplots(2, num_slices)
    colors = ['coolwarm'] * num_slices + ['coolwarm_r'] * num_slices
    for sl, ax, c in zip(slices, axes.flat, colors):
        mp, date = sl
        ax.imshow(ele, cmap='gray_r')
        ax.imshow(mp, cmap=c)
        ax.set_title(f'{date} BP')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    fig.set_size_inches(12, 6)
    fig.tight_layout()
    plt.savefig('img/slices.jpeg', dpi=300)


def plot_maps(maps, filename=None):
    fig, axes = plt.subplots(1, 2)
    for mp, ax in zip(maps, axes.flat):
        ax.imshow(ele, cmap='gray_r')
        im = ax.imshow(mp, cmap='viridis')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    plt.axis('equal')
    divider = make_axes_locatable(ax)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical',
                        fraction=0.046, pad=0.04)
    cbar.set_label('sim BP')
    fig.set_size_inches(12, 6)
    plt.savefig('img/maps.jpeg', dpi=300, bbox_inches='tight')


def test_model(params):
    time_slices = []
    maps = []
    scores = []

    num_gen = (params['start_date'] - 500) // STEP

    for i in [0, 1]:
        model = Model(**params, forest=i)
        time_slices += model.run(num_gen)
        maps.append(model.write())
        scores.append(model.score(SITES))

    return time_slices, maps, scores


def main():
    params = {'start_date': 5800,
              'start_coords': (-61.96, -10.96),
              'r': 0.025,
              'e_K': 0.25}

    time_slices, maps, scores = test_model(params)
    plot_time_slices(time_slices)
    plot_maps(maps)
    plot_graphs(scores)


if __name__ == '__main__':
    main()
