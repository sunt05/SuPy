.. _tutorial_index:

Tutorials
==================================


To familiarise users with SuPy urban climate modelling and to demonstrate the functionality of SuPy, we provide the following tutorials in `Jupyter notebooks <https://jupyter.org/>`_:

.. toctree::
  :maxdepth: 1

  quick-start
  impact-studies-parallel
  external-interaction

.. Note::
    1. The Anaconda distribution is suggested as the scientific Python 3 environment for its completeness in necessary packages. Please follow the official guide for its `installation <https://docs.anaconda.com/anaconda/install/>`__.
    2. Users with less experience in Python are suggested to go through the following section first before using SuPy.

Python 101 before SuPy
----------------------

Admittedly, this header is somewhat misleading: given the enormity of Python, it's more challenging to get this section *correct* than coding SuPy per se. As such, here a collection of data analysis oriented links to useful Python resources is provided to help novices start using Python and **then** SuPy.

- `The gist of Python <https://medium.com/@louwjlabuschagne/the-gist-of-python-ff5cc05c3318>`_: a quick introductory blog that covers Python basics for data analysis.

- Jupyter Notebook: Jupyter Notebook provides a powerful notebook-based data analysis environment that SuPy users are strongly encouraged to use. Jupyter notebooks can run in browsers (desktop, mobile) either by easy local configuration or on remote servers with pre-set environments (e.g., `Google Colaboratory <https://colab.research.google.com>`_, `Microsoft Azure Notebooks <https://notebooks.azure.com>`_). In addition, Jupyter notebooks allow great shareability by incorporating source code and detailed notes in one place, which helps users to organise their computation work.

  - Installation

    Jupyter notebooks can be installed with pip on any desktop/server system and open .ipynb notebook files locally:

      .. code-block:: shell

        python3 -m pip install jupyter -U

  - Extensions: To empower your Jupyter Notebook environment with better productivity, please check out the `Unofficial Jupyter Notebook Extensions <https://jupyter-contrib-nbextensions.readthedocs.io/en/latest/>`_. Quick introductory blogs can be found `here <https://towardsdatascience.com/jupyter-notebook-extensions-517fa69d2231>`_ and `here <https://towardsdatascience.com/bringing-the-best-out-of-jupyter-notebooks-for-data-science-f0871519ca29>`_.


- pandas: `pandas` is heavily used in SuPy and thus better understanding of pandas is essential in SuPy workflows.

  - Introductory blogs:

    * `Quick dive into Pandas for Data Science <https://towardsdatascience.com/quick-dive-into-pandas-for-data-science-cc1c1a80d9c4>`_: introduction to pandas.
    * `Basic Time Series Manipulation with Pandas <https://towardsdatascience.com/basic-time-series-manipulation-with-pandas-4432afee64ea>`_: pandas-based time series manipulation.
    * `Introduction to Data Visualization in Python <https://towardsdatascience.com/introduction-to-data-visualization-in-python-89a54c97fbed>`_: plotting using pandas and related libraries.

  - A detailed tutorial in Jupyter Notebooks:

    * `Introduction to pandas <https://github.com/fonnesbeck/Bios8366/blob/master/notebooks/Section2_1-Introduction-to-Pandas.ipynb>`_
    * `pandas fundamentals <https://github.com/fonnesbeck/Bios8366/blob/master/notebooks/Section2_2-Pandas-Fundamentals.ipynb>`_
    * `Data Wrangling with pandas <https://github.com/fonnesbeck/Bios8366/blob/master/notebooks/Section2_3-Data-Wrangling-with-Pandas.ipynb>`_
