# influenza-global-reassortment
Jupyter notebooks and data - reproducible analysis from reassortment paper.

To re-run these analyses, run the code in the following notebooks:

- `Make Figure - Known Transmission Patterns.ipynb`
- `Make Figure - Proportion Reassortant.ipynb`
- `Make Figure - PWI Histogram.ipynb`
- `Make Supplementary Data - Edge List and Node List.ipynb`

For a tutorial on how to run Jupyter notebooks, see the [Jupyter project][1] website.

Software requirements:

- Python 3.4 (minimum)
- `networkx`
- `pandas`
- `numpy`
- `seaborn`
- `matplotlib`
- `biopython`
- `joblib`
- `python-Levenshtein`

To visualize lineage traces of the 2009 H1N1 or 2013 H7N9 viruses, do the following:

1. Clone the repository to disk.
2. Navigate to the directory `supp_data/viz`
3. Start a web server in the directory. One option is to use `python -m http.server 8000`. This starts a local web server in the directory with port 8000.
4. In your browser, type: `localhost:8000/h1n1.html`. Replace the port number accordingly; to visualize H7N9, use `h7n9.html` instead.

[1]: http://jupyter.org