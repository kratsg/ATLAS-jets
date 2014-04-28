import pickle
import re
import glob
import os
import numpy as np
import matplotlib.pyplot as pl

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

#this is a wrapper around file strings to ensure the directory exists
#       could use decorators...
def write_file(f):
  ensure_dir(f)
  return f

seeds = [6, 10, 15, 20, 25, 30, 35]
#seeds = [15]
triggers = [0, 50, 100, 150, 200, 250, 300, 350, 400]
#triggers = [100]

p = re.compile(".*nsj\((.*?)\).*")

figsize = (16, 12)
labelsize = 28
titlesize = 36
legendprop = prop={'size':labelsize}
legendloc = 4


for seed,trigger in [(s,t) for s in seeds for t in triggers]:
  files = glob.glob("TurnOnCurves/Step4/pickle/matched_turnon_differential_step4_seed%d_trigger%d_nsj*" % (seed, trigger) )

  pl.figure(figsize=figsize)
  for filename in files:
    data = pickle.load(file(filename))
    pl.errorbar(data['xdata'][data['nonzero_bins']], data['ydata'], yerr=data['yerr'], label=p.match(filename).groups()[0], linewidth=2, ecolor='black')

  pl.ylim((0.0,1.2))
  pl.xlim((0.0,450.0))
  pl.title("Seed Cut %d GeV, Trigger Cut %d GeV, for Step 4" % (seed, trigger), fontsize=titlesize)
  pl.xlabel("offline $p_T^\mathrm{jet}$ [GeV]", fontsize=labelsize)
  pl.ylabel("Jet Efficiency - Differential", fontsize=labelsize)
  pl.grid(True)
  pl.legend(prop=legendprop, loc=legendloc)
  pl.savefig( write_file("TurnOnCurves/merged/merged(nsj)_seed%d_trigger%d.png" % (seed, trigger)) )
  pl.close()

  files = glob.glob("TurnOnCurves/Step4/pickle/matched_J100_curve_step4_seed%d_trigger%d_nsj*" % (seed, trigger) )

  pl.figure(figsize=figsize)
  for filename in files:
    data = pickle.load(file(filename))
    pl.errorbar(data['xdata'][data['nonzero_bins']], data['ydata'], yerr=data['yerr'], label=p.match(filename).groups()[0], linewidth=2, ecolor='black')

  pl.ylim((0.0,1.2))
  pl.xlim((0.0,450.0))
  pl.title("Seed Cut %d GeV, Trigger Cut %d GeV, for Step 4" % (seed, trigger), fontsize=titlesize)
  pl.xlabel("offline $p_T^\mathrm{jet}$ [GeV]", fontsize=labelsize)
  pl.ylabel("Event Efficiency - Differential", fontsize=labelsize)
  pl.grid(True)
  pl.legend(prop=legendprop, loc=legendloc)
  pl.savefig( write_file("TurnOnCurves/merged/merged(nsj)_J100_seed%d_trigger%d.png" % (seed, trigger)) ) 
  pl.close()

