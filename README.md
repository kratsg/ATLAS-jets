Repo containing code for analyzing jets from ATLAS

Current Features
====================
Run it simply
```
$ python
>>> from summedAreaTable import *
```

- wraps points around in \phi but not in \eta
- computes efficiency for leading jets and for gFEX module using `computeEfficiency()`
- jets are represented by 2D gaussians in eta-phi of fixed radii
- can plot example events using `showEvents($d)`

### To Do
- put it in a class/module format
- freaking optimize it

Dependencies
=====================
- ROOT 5.32, 5.34 (and PyROOT bindings)
- NumPy 1.6, 1.7, 1.8
- Python 2.6, 2.7
- GNU/Linux, Mac OS
- [root\_numpy](http://rootpy.github.io/root_numpy/install.html)
