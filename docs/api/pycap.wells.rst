Wells
-----

Two main classes are in this module: `Well` and `WellResponse` objects. `Well` objects contain all the information about a pumping well and the aquifer and stream-related properties necessary to calculate responses. Those responses are either drawdown or depletion. One or more `WellResponse` objects is constructed for each `Well`, with properties and methods necessary to calculate and keep track of specific responses.

The `AnalysisProject` object is one way to combine multiple wells and responses into a consolidated workflow. Examples show lower-level options to put together workflows.


.. automodule:: pycap.wells
   :members:
   :show-inheritance: