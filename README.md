# SuPy

[![Python Version Support Status](https://img.shields.io/pypi/pyversions/supy.svg)](https://pypi.org/project/supy)
[![Latest Version Status](https://img.shields.io/pypi/v/supy.svg)](https://pypi.org/project/supy)
[![Downloads](https://pepy.tech/badge/supy)](https://pepy.tech/project/supy)
[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sunt05/SuPy/master)

[![Build Status](https://dev.azure.com/sunt05/SuPy/_apis/build/status/sunt05.SuPy?branchName=master)](https://dev.azure.com/sunt05/SuPy/_build/latest?definitionId=11?branchName=master)
[![Documentation Status](https://readthedocs.org/projects/supy/badge/?version=latest)](https://supy.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2574404.svg)](https://doi.org/10.5281/zenodo.2574404)


[**SU**EWS](https://suews-docs.readthedocs.io) that speaks **Py**thon

## Installation

SuPy requires `python` 3.5+ and can be installed with `pip` in command line prompt:

```shell
python3 -m pip install supy --upgrade
```

## Quickstart

Once installed, `supy` can be quickly started to get [SUEWS](https://suews-docs.readthedocs.io) simulations done:

```python
import supy as sp

# load sample data
df_state_init, df_forcing = sp.load_SampleData()

# run supy/SUEWS simulation
df_output, df_state_end = sp.run_supy(df_forcing, df_state_init)

# plot results
res_plot = df_output.SUEWS.loc[1, ['QN', 'QF', 'QS', 'QE', 'QH']]
res_plot.loc['2012 6 4':'2012 6 6'].resample('30T').mean().plot()
```

The above code will produce a plot of surface energy balance components as follows:

![sample plot](https://github.com/sunt05/SuPy/raw/master/sample_plot.png)

## Tutorial

Please check out [more SuPy tutorials here!](https://supy.readthedocs.io/en/latest/tutorial/tutorial.html)
