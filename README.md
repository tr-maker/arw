# arw
Code to compute the stationary distribution for the activated random walk on a connected simple graph with one sink vertex.

``stationary_dist.py`` performs the computation and saves the results in a pickle and text files.

``stationary_dist_analysis.py`` reads in the pickle and analyzes it.

``solver.py`` and ``progressbar.py`` are helper programs used for ``stationary_dist.py`` and are created by Hannah Cairns.

The ``data`` folder contains stationary distributions and related information for some small graphs.
