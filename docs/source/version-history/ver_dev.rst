.. _new_latest:

.. _new_ver_dev:

Version ver_dev
======================================================

Compatibility Improvement.

- **New**

  1. Added version info function: `show_version`.
  2. Added command line tools:

    - `SUEWS-table-converter`: convert input tables
        from older versions to newer ones (one-way only).

    - `SUEWS`: SuPy wrapper to mimic SUEWS-binary-based simulation.


- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  1. Fixed a bug in writing out multi-grid output files
  caused by incorrect dropping of temporal information by pandas .

- **Known issue**

  None
