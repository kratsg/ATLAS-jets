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

Notes
======================
Two functions
- `showEvents($d=25)`
- `computeEfficiency()`

## `showEvents($d)` ##
`$d` is how many events to generate for (default is set at `$d=25`), but you can pass in `showEvents(5)` and it will generate files with filenames of the pattern
`<numInputJets>-jets_top-<topJetEnergy>-GeV_next-<nextJetEnergy>-GeV.png`

In the terminal, this function outputs the triggerable jets that are plotted as well as the triggered jets in a format like
```
      making plot for event
            (insert triggerable jets here)
Triggered jets
      (insert triggered jets here)
```

## `computeEfficiency()` ##
This generates two efficiency plots, one is a `jets_efficiency` which is for jets that lie outside the gFEX algorithm window. The other is a `gFEX_efficiency` which is for jets whose centroid lies only within the algorithm window.

In the terminal, this function outputs how many events it has computed efficiency for in increments of 25 events. It is rather slow at a resolution of `0.05` -- for demonstration, I would change line 51 by adding a `*4` at the end so the resolution is at `0.2x0.2` [it'll compute the efficiency fast] -- but the resolution is good at `0.05` for making images.

