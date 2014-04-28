Repo containing code for analyzing jets from ATLAS

Cloning the Repo
===================
You will want to run `git clone --recursive git@github.com:kratsg/atlas_studies` so that it initializes the submodule containing the python code for `atlas_jets`. Otherwise, just do a standard `git clone` and then run `git submodule init; git submodule update` afterwards.


Current Features
====================
See `test_atlas_jets.py` for actual code that runs. Everything is packaged up.

Dependencies
=====================
- GNU/Linux, Mac OS 10.7, 10.8, 10.9
- Python 2.6, 2.7
- ROOT 5.32+ [with PyROOT bindings](http://root.cern.ch/drupal/content/pyroot)
- [NumPy 1.6+](http://www.numpy.org/)
- [Matplotlib 1.3+](http://matplotlib.org/)
- [root\_numpy](http://rootpy.github.io/root_numpy/install.html)
