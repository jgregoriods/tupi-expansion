# Simulating the Tupi expansion in South America

Jonas Gregorio de Souza<br/>
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7879--4531-green.svg)](https://orcid.org/0000-0001-6032-4443)<br/>

<p>This repository contains the R and python code for simulating demic diffusion of tropical forest farmers in lowland South America, namely the Tupi language family, and the role of climate change/forest expansion in shaping their dispersal.</p>

<p>For details, please see:</p>
<p>Souza J.G., Noelli F.S., Madella M. <i>Reassessing the role of climate change in the Tupi expansion.</i></p>

## Usage

<p>The python model implements a cellular automaton simulating population growth and dispersal in a dynamic environment. Every run starts with a single settled cell, whose population is initialised at K, immediately dispersing to the neighbouring cells. Each time step, for every cell that is inhabited, population growth and dispersal methods are applied, and simulated dates of arrival are recorded. Models can be run in which all land cells may be settled, and others in which settlement is restricted to tropical moist forest.</p>

<p>To run:</p>

<pre><code>>>> python3 python/main.py
</pre></code>

<p>This will run a number of models with different parameters and calculate the scores based on the difference from the empirical 14C dates. The scores will be written to <code>results/ca_scores.csv</code> and the simulated dates with the default parameters will be written to <code>results/sim_dates.csv</code>. Time slices (every 30 generations) and the final simulated arrival times will be written as rasters to the folder <code>results/rasters</code>.</p>

<p>The R code implements an equation-based model of front propagation with an underlying cost surface. The cost of traversing different biomes can be adjusted for the cost surface, e.g. slowing down the expansion in non-forested biomes. All simulations start 5000 BP at the coordinates of the Tupi site with the earliest date.</p>

<p>To run, start R and change to the respective directory:</p>

<pre><code>> setwd("R")
> source("main.R")
</pre></code>

<p>This will run a number of models with different parameters and values for the cost surface. Scores will be written to <code>results/ebm_scores.csv</code>. The figures used in the article (arrival times under different models, scatterplots of real and simulated dates vs distances from origin, time slices) will also be written to folder <code>img</code>.</p>