Contributing to HiCap Analysis
==============================

This is an open-source project, so contributions
are welcome.  The proposed approach is to clone
repository, work in a branch, test the feature
you are adding and then merge into master.

Try to add tests to the test directory.  Make sure they
run locally by running (from a shell in the main directory)

::

    coverage run -m pytest
    coverage report -m 

Use docstrings.  Documentation can be added to the
source directory in `rst` files.  To have sphinx
pull in docstrings, add a `rst` file for the
module to the source/api directory 
and an entry to `index.rst` in 
the docs/source/api directory.

In `index.rst`

::

    .. toctree::
       pycap.name_of_module-no_py_extension


The module `rst` file looks like

::

    Module name
    -----------

    Long description if desired....  the automodule commands pull in the docstrings.

    .. automodule:: pycap.utilities
          :members:
          :show-inheritance:


:: 