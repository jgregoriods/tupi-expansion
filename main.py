import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import seaborn as sns

from itertools import product
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from model import *


sns.set(font_scale=1.5)

SITES = pd.read_csv('sites/tupi_filtered.csv')
ele = np.loadtxt('layers/ele.asc', skiprows=6)
ele[ele < 1] = np.nan
ele[ele >= 1] = 1
sam = gpd.read_file('shp/south_america.shp')

blue = '#f9b641ff'
red = '#cc6a70ff'
green = '#6b4596ff'
blue_cmap = ListedColormap([blue])
red_cmap = ListedColormap([red])
green_cmap = ListedColormap([green])
white = ListedColormap(['white'])
black = ListedColormap(['black'])


class ModelTest:
    def __init__(self, start_dates, start_coords, rs, e_Ks):
        self.start_dates = start_dates
        self.start_coords = start_coords
        self.rs = rs
        self.e_Ks = e_Ks

        self.sites = SITES

        self.models = []

    def run_models(self):
        for i in list(product(self.start_dates, self.start_coords, self.rs, self.e_Ks)):
            params = {'start_date': i[0],
                      'start_coord': i[1],
                      'r': i[2],
                      'e_K': i[3]}
            num_gen = (i[0] - 500) // STEP
            for j in ['null', 'dry', 'moist']:
                model = Model(**params, forest=j)
                model.run(num_gen)
                model.get_score(SITES)

                self.models.append(model)
                self.sites[f'sim_{j}'] = model.sites['sim_dates']

    def plot_maps(self, filepath=None):
        maps = [model.arrival_times for model in self.models]
        fig, axes = plt.subplots(1, len(maps))
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
        fig.set_size_inches(18, 6)
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def plot_time_slices(self, filepath=None):
        slices = [i for model in self.models for i in model.slices]
        num_slices = len(slices) // 3
        fig, axes = plt.subplots(3, num_slices)
        #colors = [blue_cmap] * 6 + [red_cmap] * 6 + [green_cmap] * 6
        #for sl, ax, c in zip(slices, axes.flat, colors):
        for sl, ax in zip(slices, axes.flat):
            mp, date = sl
            ax.imshow(ele, cmap=white, extent=[-81.34,-34.79,-55.92,12.47])
            ax.imshow(mp, cmap=black, extent=[-81.34,-34.79,-55.92,12.47], zorder=3, vmax=2)
            ax.set_title(f'{date} BP')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
        fig.set_size_inches(12, 12)
        fig.tight_layout()
        if filepath is not None:
            plt.savefig(filepath, dpi=300)
        else:
            plt.show()

    def plot_graphs(self, filepath=None):
        scores = [model.sites for model in self.models]
        fig, axes = plt.subplots(1, len(scores))
        #for score, ax, c in zip(scores, axes.flat, [blue, red, green]):
        for score, ax in zip(scores, axes.flat):
            ax.scatter(SITES['dist'], SITES['bp'], c='darkgray')
            ax.scatter(score['dist'][score['sim_dates'] != 0],
                       score['sim_dates'][score['sim_dates'] != 0],
                       c='black')
            ax.set_xlabel('distance (km)')
            ax.set_ylabel('age (cal BP)')
        fig.set_size_inches(18, 6)
        fig.tight_layout()
        if filepath is not None:
            plt.savefig('img/plots.jpeg', dpi=300)
        else:
            plt.show()

    def write_sim_dates(self, filepath):
        np.savetxt(filepath, self.sites, fmt='%5s', delimiter=',',
                   header=','.join(self.sites.columns), comments='')


def main():
    mt = ModelTest([4000], [(-61.96, -10.96)], [0.028], [0.25])
    mt.run_models()
    mt.plot_maps('img/maps.jpeg')
    mt.plot_time_slices('img/time_slices.jpeg')
    mt.plot_graphs('img/graphs.jpeg')
    mt.write_sim_dates('img/sim_dates.csv')


if __name__ == '__main__':
    main()
