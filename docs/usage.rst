=====
Usage
=====
Overview
--------
The `pycap-dss` software provides a suite of analytical solutions to quantify stream depletion and drawdown due to pumping a groundwater well. In the sections below are information about individual solutions. These can be exposed and evaluated directly. A list of available functions for depletion and drawdown, respectively, can be listed from the code as follows:

.. code-block:: python

    import pycap
    # list the depletion functions
    pycap.ALL_DEP_METHODS
    # list the drawdown functions
    pycap.ALL_DD_METHODS


The code is designed in a modular way such that additional methods can be implemented and then accessed by the rest of the code. The following sections expand on these capabilities and build on one another in terms of complexity. The :ref:`Individual Solutions` section describes how each analytical solution can be accessed directly and executed. The :ref:`Well Object` section describes the `Well` class that is a wrapper around the depletion and drawdown solution and also can be initialized with general properties about wells and responses (e.g. streams or drawdown locations). A single instance of a `Well` object represents a pumping well and can include aquifer properties, pumping schedules, and various responses at various locations. Finally, the :ref:`Analysis Project` section describes an example project as implemented by the Wisconsin Department of Natural Resources in which multiple existing and/or proposed pumping wells with multiple responses can be evaluated in a coordinated single project. This represents one example, but not the only way, to coordinate a potentially complex area with many wells and responses and to systematically report the results.

Individual Solutions
--------------------
The: `Solutions_Demonstration <examples/Solutions_Demonstration.ipynb>` notebook demonstrated the use of individual analytical solutions in `pycap-dss`.


Well Object
-----------

    - allows for superposition in time making time series pumping
    - multiple responses
    - changing the methods and apportionment

Analysis Project
----------------

    - an example wrapping it all up together




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

    
