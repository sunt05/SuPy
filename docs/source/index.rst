.. _index_page:

SuPy: SUEWS that speaks Python
----------------------------------------------------

.. image:: https://img.shields.io/pypi/pyversions/supy.svg
    :target: https://pypi.org/project/supy
    :alt: Python Version Support Status

.. image:: https://img.shields.io/pypi/v/supy.svg
    :target: https://pypi.org/project/supy
    :alt: Latest Version Status

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/sunt05/SuPy/master
    :alt: Binder Status


.. image:: https://readthedocs.org/projects/supy/badge/?version=latest
    :target: https://supy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://dev.azure.com/sunt05/SuPy/_apis/build/status/sunt05.SuPy?branchName=master
    :target: https://dev.azure.com/sunt05/SuPy/_build/latest?definitionId=11?branchName=master
    :alt: Build Status


.. caution::

    This site is under construction. All information might NOT be accurate and are subject to rapid change.



- **What is SuPy?**

    SuPy is a Python-enhanced urban climate model
    with `SUEWS`_ as its computation core.

    The scientific rigour in SuPy results is thus gurranteed by SUEWS
    (see :ref:`SUEWS publications <Recent_publications>` and
    :ref:`Parameterisations and sub-models within SUEWS`).

    Meanwhile, the data analysis ability of SuPy is greatly enhanced
    by `the Python-based SciPy Stack <https://scipy.org>`_, notably `numpy`_ and `pandas`_.


.. _SUEWS: https://suews-docs.readthedocs.io/en/latest/
.. _numpy: https://www.numpy.org
.. _pandas: http://pandas.pydata.org/


- **How to get SuPy?**

  SuPy is available on all major platforms (macOS, Windows, Linux) for Python 3.5+
  via `PyPI <https://pypi.org/project/supy/>`_:

  .. code-block:: shell

    python3 -m pip install supy --upgrade

- **How to use SuPy?**

    * Please follow :ref:`Quickstart of SuPy`
    and :ref:`other tutorials <tutorial_index>`.

    * Please see `api` for details.

- **How to contribute to SuPy?**

    * Add your development via `Pull Request <https://help.github.com/articles/about-pull-requests/>`_
    * Report issues via the `GitHub page <https://github.com/sunt05/SuPy/issues/new?template=issue-report.md>`_.
    * Provide suggestions and feedback.

.. toctree::
  :hidden:
  :maxdepth: 1

  tutorial/tutorial
  data-structure/supy-io
  api
  version-history/version-history

