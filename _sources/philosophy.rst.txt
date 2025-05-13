Philosophy
==========

Motivation
----------
Analytical solutions for streamflow depletion or 
drawdown from pumping wells are being used to quickly screen
proposed high-capacity pumping impacts.  Python implementations
of these solutions can be easier for analysts to access compared
to FORTRAN-based solutions.  Python also can be faster and provide
a more robust option compared to spreadsheet implementations 
of the solutions.


What pycap does
-------------------------------------------
The package defines classes for analysis of the potential 
responses of water levels and stream depletion resulting from
pumping of a high-capacity well. The user performs analysis for a 
specific location for analysis of streamflow depletion or
drawdown.  The project class also allows for the analysis
of several pumping wells at the same time.

The underlying analytical solutions also are available
as internal functions of the wells class using the 
underscore notation to indicate that these are typically
called from within the class.  For example direct access
to the Theis (1935) drawdown solution is provided by
`wells.theis(T,S,time,dist,Q)` function.


What pycap doesn't do
--------------------------------------------------------
The modules assembled here will be used in notebooks or scripts
to actually help build the final model.