Repo containing code for analyzing jets from ATLAS

Current Features
====================

- adds jets to a grid that satisfy jetPt > 200. GeV and are within the grid's boundaries of `[phi, eta] = [(-3.2,3.2), (0,3.2)]`
- reports metrics about centroid coordinates that did not get added to grid
- reports how many jets are in an event
- reports how many jets for an event did not trigger (above threshold of 200 GeV)
- reports how many jets for an event added to grid
- reports how many points for a jet are added to grid
- wraps points around in \phi but not in \eta
- computes efficiency for leading jets and for gFEX module
- gaussians and circular mesh grids are added

### To Do
- put it in a class/module format

Dependencies
=====================
- ROOT 5.32, 5.34
- NumPy 1.6, 1.7, 1.8
- Python 2.6, 2.7
- GNU/Linux, Mac OS
- root\_numpy [root\_numpy](http://rootpy.github.io/root_numpy/install.html)
