import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import seaborn as sns
sns.set(font_scale=1.5)

from itertools import product
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from model import *


SITES = pd.read_csv('sites/tupi_filtered_2.csv')
ele = np.loadtxt('layers/ele.asc', skiprows=6)
ele[ele < 1] = np.nan
ele[ele >= 1] = 1
sam = gpd.read_file('shp/south_america.shp')

blue = '#f9b641ff'
red = '#cc6a70ff'
blue_cmap = ListedColormap([blue])
red_cmap = ListedColormap([red])
white = ListedColormap(['white'])


def plot_graphs(scores, filepath=None):
    fig, axes = plt.subplots(1, 2)
    cs = ['#0449e1', '#c00000']
    for score, ax, c in zip(scores, axes.flat, [blue, red]):
        ax.scatter(SITES['dist'], SITES['bp'], c='black')
        ax.scatter(score['dist'][score['sim_dates'] != 0],
                   score['sim_dates'][score['sim_dates'] != 0],
                   c=c)
        ax.set_xlabel('distance (km)')
        ax.set_ylabel('age (cal BP)')
    fig.set_size_inches(12, 6)
    fig.tight_layout()
    if filepath is not None:
        plt.savefig('img/plots.jpeg', dpi=300)
    else:
        plt.show()


def plot_time_slices(slices, filepath=None):
    num_slices = len(slices) // 2
    fig, axes = plt.subplots(2, num_slices)
    colors = [blue_cmap] * 6 + [red_cmap] * 6
    for sl, ax, c in zip(slices, axes.flat, colors):
        mp, date = sl
        ax.imshow(ele, cmap=white, extent=[-81.34,-34.79,-55.92,12.47])
        ax.imshow(mp, cmap=c, extent=[-81.34,-34.79,-55.92,12.47], zorder=3, vmax=2)
        ax.set_title(f'{date} BP')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    fig.set_size_inches(12, 6)
    fig.tight_layout()
    if filepath is not None:
        plt.savefig(filepath, dpi=300)
    else:
        plt.show()


def plot_maps(maps, filepath=None):
    fig, axes = plt.subplots(1, 2)
    for mp, ax in zip(maps, axes.flat):
        mp[mp==0] = np.nan
        ax.imshow(ele, cmap=white, extent=[-81.34,-34.79,-55.92,12.47])
        im = ax.imshow(mp, cmap='viridis', extent=[-81.34,-34.79,-55.92,12.47], zorder=3)
        #ax.xaxis.set_ticklabels([])
        #ax.yaxis.set_ticklabels([])
    plt.axis('equal')
    divider = make_axes_locatable(ax)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical',
                        fraction=0.046, pad=0.04)
    cbar.set_label('sim BP')
    fig.set_size_inches(12, 6)
    if filepath is not None:
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
    else:
        plt.show()


def test_model(params):
    time_slices = []
    maps = []
    scores = []
    rmse = SITES

    num_gen = (params['start_date'] - 500) // STEP

    for i in [0, 1]:
        model = Model(**params, forest=i)
        model.run(num_gen)
        model.get_score(SITES)

        time_slices += model.slices
        maps.append(model.arrival_times)
        scores.append(model.sites)
        rmse[f'sim{i}'] = model.sites['sim_dates']

    return time_slices, maps, scores, rmse


def main():
    params = {'start_date': 5800,
              'start_coords': (-61.96, -10.96),
              'r': 0.02,
              'e_K': 0.7}

    time_slices, maps, scores, rmse = test_model(params)

    plot_time_slices(time_slices, 'img/time_slices.jpeg')
    plot_maps(maps, 'img/maps.jpeg')
    plot_graphs(scores, 'img/graphs.jpeg')
    np.savetxt('img/rmse.csv', rmse, fmt='%5s', delimiter=',',
               header=','.join(rmse.columns), comments='')


if __name__ == '__main__':
    main()
