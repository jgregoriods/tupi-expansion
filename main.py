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

        self.sites_df = pd.DataFrame(columns=['id', 'forest', 'r', 'e_K', 'dist', 'sim_dates'])

        self.models = []

    def run_models(self):
        model_id = 1
        for i in list(product(self.start_dates, self.start_coords, self.rs, self.e_Ks)):
            params = {'start_date': i[0],
                      'start_coord': i[1],
                      'r': i[2],
                      'e_K': i[3]}
            num_gen = (i[0] - 500) // STEP
            for j in ['null', 'moist']:
                model = Model(**params, forest=j)
                model.run(num_gen)
                model.get_score(SITES)

                self.models.append(model)

                df = pd.DataFrame(columns=['id', 'forest', 'r', 'e_K', 'dist', 'sim_dates'])
                df['dist'] = model.sites['dist']
                df['sim_dates'] = model.sites['sim_dates']
                df['id'] = model_id
                df['forest'] = model.forest
                df['r'] = model.r
                df['e_K'] = model.e_K

                self.sites_df = pd.concat([self.sites_df, df])

            model_id += 1

    def plot_maps(self, filepath=None):
        extent = [-81.34,-34.79,-55.92,12.47]
        maps = [model.arrival_times for model in self.models]
        fig, axes = plt.subplots(len(maps) // 2, 2)
        cmap = plt.get_cmap('viridis', 10)
        for mp, ax in zip(maps, axes.flat):
            mp[mp==0] = np.nan
            ax.imshow(ele, cmap=white, extent=extent)
            im = ax.imshow(mp, cmap=cmap, extent=extent, zorder=3)
        plt.axis('equal')
        divider = make_axes_locatable(ax)
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical', fraction=0.046, pad=0.04)
        cbar.set_label('sim BP')
        fig.set_size_inches(12, len(maps)*3)
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def plot_time_slices(self, filepath=None):
        if len(self.models) > 2:
            return
        extent = [-81.34,-34.79,-55.92,12.47]
        slices = [i for model in self.models for i in model.slices]
        num_slices = len(slices) // 2
        fig, axes = plt.subplots(2, num_slices)
        for sl, ax in zip(slices, axes.flat):
            mp, date = sl
            ax.imshow(ele, cmap=white, extent=extent)
            ax.imshow(mp, cmap=black, extent=extent, zorder=3, vmax=2)
            ax.set_title(f'{date} BP')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
        fig.set_size_inches(12, 6)
        fig.tight_layout()
        if filepath is not None:
            plt.savefig(filepath, dpi=300)
        else:
            plt.show()

    def plot_graphs(self, filepath=None):
        fig, axes = plt.subplots(1, 2)
        models = ['null', 'moist']
        for i, ax, forest in zip(range(1, len(self.models) + 1), axes.flat, models):
            sites_to_plot = self.sites_df.loc[(self.sites_df['forest'] == forest) &
                                              (self.sites_df['sim_dates'] != 0)]
            ax.scatter(SITES['dist'], SITES['bp'], c='darkgray')
            if len(self.models) > 2:
                ax.scatter(sites_to_plot['dist'], sites_to_plot['sim_dates'],
                           c=sites_to_plot['r'], cmap='viridis',
                           s=sites_to_plot['e_K'] * 100)
            else:
                ax.scatter(sites_to_plot['dist'], sites_to_plot['sim_dates'], c='black')
            ax.set_xlabel('distance (km)')
            ax.set_ylabel('age (cal BP)')
        fig.set_size_inches(12, 6)
        fig.tight_layout()
        if filepath is not None:
            plt.savefig(filepath, dpi=300)
        else:
            plt.show()

    def write_sim_dates(self, filepath):
        np.savetxt(filepath, self.sites_df, fmt='%5s', delimiter=',',
                   header=','.join(self.sites_df.columns), comments='')


def main():
    mt = ModelTest(start_dates=[5000],
                   start_coords=[(-61.96, -10.96)],
                   rs=[0.025],
                   e_Ks=[0.25])
    mt.run_models()
    mt.plot_maps('img/maps.pdf')
    mt.plot_time_slices('img/time_slices.pdf')
    mt.plot_graphs('img/graphs.pdf')
    mt.write_sim_dates('img/sim_dates.csv')


if __name__ == '__main__':
    main()
