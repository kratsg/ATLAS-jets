from atlas_jets import *
import ROOT
import numpy as np
import matplotlib.pyplot as pl
import root_numpy as rnp
import pickle

# define helper functions - also a source of parallelization!
def compute_jetDistance(jet1, jet2):
  return ((jet1.eta - jet2.eta)**2. + (jet1.phi - jet2.phi)**2.)**0.5

def match_jets(oJets=[], tJets=[]):
  # dR as a change in distance
  dR = 1.0
  if len(tJets) == 0:
    return np.array([[oJet,gTowers.Jet()] for oJet in oJets])
  # we want to match the closest gTower jet for every offline jet
  matched_jets = []
  for oJet in oJets:
    distances = np.array(map(lambda tJet: compute_jetDistance(tJet, oJet), tJets))
    energies = np.array(map(lambda tJet: tJet.E/np.cosh(tJet.eta), tJets))
    # return jet with highest ET within dR
    if np.where(distances <= dR)[0].size == 0:
      closest_tJet = gTowers.Jet()
    else:
      max_energy_in_distance = np.amax(energies[np.where(distances <= 1.0)])
      index_jet = np.where(energies == max_energy_in_distance)[0][0]
      closest_tJet = tJets[index_jet]
    matched_jets.append([oJet,closest_tJet])
  return matched_jets

#no 22 or 23
def run_code(page_num = 0):

  filename = '/Users/kratsg/Desktop/LArStudies/PileupSkim_TTbar_14TeV_MU80_%d.root' % page_num
  #directory = 'TTbar_14TeV_MU80'
  tree = 'mytree'

  rootfile = ROOT.TFile.Open(filename)

  #set total number of events
  total_num_events = int(rootfile.Get('%s' % (tree)).GetEntries())

  # seed_ETthresh = 15
  seed_ETthresh = 15.

  #set seed cuts
  seed_filter = gTowers.SeedFilter(ETthresh = seed_ETthresh, numSeeds = 100)
  #column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
  offline_column_names = ['jet_AntiKt10LCTopoTrimmedPtFrac5SmallR30_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi', 'TrimmedSubjetsPtFrac5SmallR30_nsj', 'Tau1', 'Tau2', 'Tau3', 'SPLIT12', 'SPLIT23', 'SPLIT34']]
  gTower_column_names = ['gTower%s' % col for col in ['E', 'Et', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]

  #paired jets initialization
  paired_jets = []

  # main loop that goes over the file
  for event_num in range(total_num_events):
    if event_num % 100 == 0:
      print "doing event_num=%d for page_num: %d" % (event_num, page_num)
    # pull in data row by row
    data = rnp.root2rec(filename, treename='%s' % (tree), branches=offline_column_names + gTower_column_names, start=(event_num), stop=(event_num+1))
    print '\t%d: loaded data' % page_num
    oEvent = OfflineJets.Event(event=[data[col][0] for col in offline_column_names])

    # if there are no offline jets, we skip it
    if len(oEvent.jets) == 0:
      continue

    '''can use seed_filter on an event by event basis'''
    # max number of seeds based on number of offline jets
    #seed_filter = gTowers.SeedFilter(numSeeds = len(oEvent.jets))
    print '\t%d: building tEvent' % page_num
    tEvent = gTowers.TowerEvent(event=[data[col][0] for col in gTower_column_names], seed_filter = seed_filter)
    print '\t%d: done' % page_num
    # build up the first two histograms using just the gTower data
    # note, we have np.histogram(...)[0] since we only need the hist data

    tEvent.get_event()
    #paired_jets = match_jets(oJets=oEvent.jets, tJets=tEvent.filter_towers())
    paired_jets.append( match_jets(oJets=oEvent.jets, tJets=tEvent.event.jets) )

  '''at this point, we've processed all the data and we just need to dump it'''
  filename_ending = 'seed%d_unweighted_page%d' % (seed_ETthresh, page_num)
  pickle.dump(paired_jets, file('matched_jets_%s.pkl' % filename_ending, 'w+') )
  return len(paired_jets)

class Copier(object):
  def __init__(self):
    pass
  def __call__(self, args):
    page = args
    run_code(page)

#pages = [0,1,10,11]
pages = [12,14,15,16]
#pages = [0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,23]
jobs = [(c) for c in pages]
print jobs

import multiprocessing
p = multiprocessing.Pool(processes=4)
p.map(Copier(), jobs)
