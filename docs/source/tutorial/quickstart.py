#%% [markdown]
# # Quickstart of SuPy
# 
# This quickstart demonstrates the essential and simplest workflow of `supy` in SUEWS simulation:
# 
# 1. [load SUEWS input files](#Load-SUEWS-input-files)
# 2. [run SUEWS simulation](#Run-SUEWS-simulations)
# 3. [examine SUEWS results](#Examine-SUEWS-results)
# 
# More advanced use of `supy` are available in the [tutorials](tutorial)
# 
# Before start, we need to load the following necessary packages.

#%%
import supy as sp
import pandas as pd
import numpy as np
from pathlib import Path
get_ipython().run_line_magic('matplotlib', 'inline')

#%% [markdown]
# ## Load SUEWS input files
#%% [markdown]
# First, a path to SUEWS `RunControl.nml` should be specified, which will direct `supy` to locate input files.

#%%
path_runcontrol=Path('../sample_run')/'RunControl.nml'

#%% [markdown]
# ### Model configuration
# 
# We call `sp.init_supy` to initialise a SUEWS simulation and get two `pandas` objects (note: the following names CAN be customised and are NOT fixed to the examples shown here):
# 
#  `df_state_init`: a `DataFrame` for grid-specific settings
# 
# Once loaded in, these objects CAN be modified and reused for conducting simulations that differ from the one configured via input files under the above `dir_start`.

#%%
df_state_init = sp.init_supy(path_runcontrol)

#%% [markdown]
# A sample `df_state_init` looks below (note that `.T` is used here to a nicer tableform view):

#%%
df_state_init.filter(like='method').T

#%% [markdown]
# ### Meteorological forcing
# 
# Following the convention of SUEWS, `supy` loads meteorological forcing (met-forcing) files at the grid level.
# 
# <div class="alert alert-info">
# 
# **Note:** 
#     
#     If `multiplemetfiles = 0` (i.e., all grids use the same met-forcing file) is set in `ser_mod_cfg`, the `grid` argument takes NO effect and is ignored by `supy`.
# 
# </div>

#%%
grid = df_state_init.index[0]
df_forcing = sp.load_forcing_grid(path_runcontrol, grid)

#%% [markdown]
# ## Run SUEWS simulations
#%% [markdown]
# Once met-forcing (via `df_forcing`) and initial conditions (via `df_state_init`) are loaded in, we call `sp.run_supy` to conduct a SUEWS simulation, which will return two `pandas` `DataFrame`s: `df_output` and `df_state`.

#%%
df_output, df_state = sp.run_supy(df_forcing, df_state_init)

#%% [markdown]
# ### `df_output`
# 
# `df_output` is an ensemble output collection of major SUEWS output groups, including:
#     
#     * SUEWS: the essential SUEWS output variables
#     * DailyState: variables of daily state information
#     * snow: snow output variables (effective when `snowuse = 1` set in `ser_mod_cfg`)
#     * ESTM: ESTM output variables (not implemented yet)
# 

#%%
df_output.columns.levels[0]

#%% [markdown]
# ### `df_state`
# 
# `df_state` is a `DataFrame` for holding:
#     
#    1. all model states if `save_state` is set to `True` when calling `sp.run_supy` and `supy` may run significantly slower for a large simulation;
#    2. or, only the final state if `save_state` is set to `False` (the default setting) in which mode `supy` has a similar performance as the standalone compiled SUEWS executable.
# 
# Entries in `df_state` have the same data structure as `df_state_init` and can thus be used for other SUEWS simulations staring at the timestamp as in `df_state`.

#%%
df_state.T

#%% [markdown]
# ## Examine SUEWS results
#%% [markdown]
# Thanks to the functionality inherited from `pandas` and other packages under the [PyData](https://pydata.org) stack, compared with the standard SUEWS simulation workflow, `supy` enables more convenient examination of SUEWS results by statistics calculation, resampling, plotting (and many more).
#%% [markdown]
# ### Ouptut structure
# 
# `df_output` is organised with `MultiIndex` `(grid,timestamp)` and `(group,varaible)` as `index` and `columns`, respectively.

#%%
df_output.head()

#%% [markdown]
# Here we demonstrate several typical scenarios for SUEWS results examination.
# 
# The essential `SUEWS` output collection is extracted as a separate variable for easier processing in the following sections. More [advanced slicing techniques](http://pandas.pydata.org/pandas-docs/stable/advanced.html#multiindex-advanced-indexing) are available in `pandas` documentation.

#%%
df_output_suews=df_output['SUEWS']

#%% [markdown]
# ### Statistics Calculation
# 
# We can use `.describe()` method for a quick overview of the key surface energy balance budgets.

#%%
df_output_suews.loc[:,['QN','QS','QH','QE','QF']].describe()

#%% [markdown]
# ### Plotting
# 
# Plotting is very straightforward via the `.plot` method bounded with `pandas` objects.

#%%
df_output_suews.loc[1].loc[
    '2012 6 1':'2012 6 7',
    ['QN','QS','QE','QH','QF']
].plot()

#%% [markdown]
# ### Resampling
# 
# The suggested runtime/simulation frequency of SUEWS is `300 s`, which usually results a large output and may be over-weighted for storage. To slim down the output size, we can `resample` the default output. 

#%%
df_output_suews_rsmp=df_output_suews.loc[1].resample('1h').mean()


#%%
df_output_suews_rsmp.loc[
    '2012 6 1':'2012 6 7',
    ['QN','QS','QE','QH','QF']].plot()

#%% [markdown]
# The resampled output can be outputed for a smaller file.

#%%
df_output_suews_rsmp.to_csv('suews_1h.txt',
                            sep='\t',
                            float_format='%8.2f',
                            na_rep=-999)

#%% [markdown]
# For a justified format, we use the `to_string` for better format controlling and write the formatted string out to a file.

#%%
str_out=df_output_suews_rsmp.to_string(
    float_format='%8.2f',
    na_rep='-999',
    justify='right')
with open('suews_sample.txt','w') as file_out:
    print(str_out,file=file_out)


