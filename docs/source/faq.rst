.. _faq:



FAQ
===

.. contents:: Contents
   :local:
   :backlinks: none

I cannot install SuPy following the docs, what is wrong there?
----------------------------------------------------------------

please check if your environment meets the following requirements:

1. Operating system (OS):

  a. is it 64 bit? only 64 bit systems are supported.

  b. is your OS up to date? only recent desktop systems are supported:

    - Windows 10 and above
    - macOS 10.13 and above
    - Linux: no restriction;
      If SuPy cannot run on your specific Linux distribution,
      please report it to us.

    You can get the OS information with the following code:

    .. code-block:: python

        import platform
        platform.platform()

2. Python interpreter:

  a. is your Python interpreter 64 bit?

    Check running mode with the following code:

    .. code-block:: python

        import struct
        struct.calcsize('P')*8

  b. is your Python version above 3.5?

    Check version info with the following code:

    .. code-block:: python

        import sys
        sys.version

If your environment doesn't meet the requirement by SuPy,
please use a proper environment;
otherwise, `please report your issue`__.

__ new_issue_

How do I know which version of SuPy I am using?
-----------------------------------------------

Use the following code:

.. code-block:: python

    import supy
    supy.show_version()

.. note:: `show_version` is introduced after supy v2019.5.28; previous versions doesn't have this function.



A `kernel may have died` exception happened, where did I go wrong?
------------------------------------------------------------------

The issue is highly likely due to invalid input to SuPy and SUEWS kernel.
We are trying to avoid such exceptions,
but unfortunately they might happen in some edge cases.

Please `report such issues to us`__ with your input files for debugging.
Thanks!

__ new_issue_


How can I upgrade SuPy to an up-to-date version?
------------------------------------------------
Run the following code in your terminal:

.. code-block:: python

    python3 -m pip install supy --upgrade



