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


What pycap-dss does
-------------------------------------------
The package defines classes for analysis of the potential 
responses of water levels and stream depletion resulting from
pumping of a high-capacity well. The user performs analysis for a 
specific location for analysis of streamflow depletion or
drawdown.  The project class also allows for the analysis
of several pumping wells at the same time.

The code is designed in a modular way such that new depletion or drawdown solutions can be cleanly added and then called by the higher-level objects.

The underlying analytical solutions also are available
as internal functions of the wells class using the 
underscore notation to indicate that these are typically
called from within the class.  For example direct access
to the Theis (1935) drawdown solution is provided by
`wells.theis_drawdown(T,S,time,dist,Q)` function.


What pycap-dss doesn't do
--------------------------------------------------------
This package does not include all possible analytical solutions for drawdown or stream depletion due to wells. This package is not a replacement for numerical or other modeling approaches in situations with complex hydrogeologic conditions.