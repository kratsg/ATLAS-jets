import matplotlib.pyplot as pl
import numpy as np
import pickle
from glob import glob
import re

def tryint(s):
  try:
    return int(s)
  except:
    return s

def alphanum_key(s):
  """ Turn a string into a list of string and number chunks.
      "z23a" -> ["z", 23, "a"]
  """
  return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
  """ Sort the given list in the way that humans expect.
  """
  l.sort(key=alphanum_key)

histogram_data = []
seedcuts = []

p = re.compile('^.*?_seed(\d+)_.*?\.pkl$')
files = glob("events_histogram_leading_trigger_jets_offline0_gTower0_seed*.pkl")
sort_nicely(files)
for filename in files:
  data = pickle.load(file(filename))
  seed_filter_ETthresh = p.match(filename).groups()
  # load data
  histogram_data.append(data)
  seedcuts.append(p.match(filename).groups()[0])
  pl.close()

pl.close()

colors = ['b','g','r','c','m','y','k','w']

pl.figure(figsize=(25.0, 15.0))
for hist,seedcut,color in zip(histogram_data, seedcuts, colors):
  values = hist['bins']
  weights = hist['values']
  #bins = hist['bins']
  bins = np.arange(0.0,4000.0,10.0)#hist[:,0]
  pl.hist(values[:-1], bins=bins, weights=weights, histtype='bar', label='%s GeV' % seedcut, alpha=0.9, color=color)

filename_ending = "-".join(seedcuts)

pl.xlabel('$E_T^{\mathrm{jet}}$ [GeV]')
pl.ylabel('Number of leading trigger jets')
pl.title('$p_T^{\mathrm{offline jet}}$ > 0 GeV, 10000 events, $E_T^{\mathrm{tower jet}}$ > 0 GeV')
pl.legend()
#pl.xlim(xlim)
pl.xlim((50.0,550.0))
pl.grid(True, which='both', linewidth=3, linestyle='--')
#pl.savefig('events_histogram_leading_trigger_jets_offline0_gTower0_seed%s_unweighted.png' % filename_ending, dpi=200)
pl.savefig('events_histogram_leading_trigger_jets_offline0_gTower0_seed%s_unweighted_CROPPED.png' % filename_ending, dpi=200)
pl.close()
