=====
Usage
=====

# Individual Solutions

# Well object
    - allows for superposition in time making time series pumping
    - multiple responses
    - changing the methods and apportionment

# Analysis Project
    - an example wrapping it all up together

# maybe geoprocessing?



Start by importing `pycap`

.. code-block:: python

    import pycap


Or import a well functions

.. code-block:: python

    import pycap.wells as wo

    ddwn = wo.theis()

Or import a specific function.

.. code-block:: python

    from pycap.wells import glover

    depl = glover()

    
