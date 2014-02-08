from atlas_jets import *
import ROOT
import numpy as np
import matplotlib.pyplot as pl
import root_numpy as rnp
import pickle
filename = '/Users/kratsg/Desktop/PileupSkim_TTbar_14TeV_MU80_10000.root' 
directory = 'TTbar_14TeV_MU80'
tree = 'mytree'

rootfile = ROOT.TFile(filename)
#set total number of events
total_num_events = int(rootfile.Get('%s/%s' % (directory, tree)).GetEntries())

#set jet cuts
offline_jetpT_threshold = 0. #[GeV]
gTower_jetET_threshold  = 0.

# define helper functions - also a source of parallelization!
def compute_jetDistance(jet1, jet2):
  return ((jet1.eta - jet2.eta)**2. + (jet1.phi - jet2.phi)**2.)**0.5

def match_jets(oJets=[], tJets=[]):
  if len(tJets) == 0:
    return np.array([[oJet,gTowers.Jet()] for oJet in oJets])
  # we want to match the closest gTower jet for every offline jet
  matched_jets = []
  for oJet in oJets:
    distances = np.array(map(lambda tJet: compute_jetDistance(tJet, oJet), tJets))
    index_closest = np.argmin(distances)
    if distances[index_closest] > 1.0:
      closest_tJet = gTowers.Jet()
    else:
      closest_tJet = tJets[np.argmin(distances)]
    matched_jets.append([oJet,closest_tJet])
  return np.array(matched_jets)

def run_code(offline_jetpT_threshold = 0., gTower_jetET_threshold = 0., seed_ETthresh = 0.):
  #set seed cuts
  seed_filter = gTowers.SeedFilter(ETthresh = seed_ETthresh, numSeeds = 1.0e5)

  #column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
  offline_column_names = ['jet_AntiKt4LCTopo_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi']]
  gTower_column_names = ['gTower%s' % col for col in ['E', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]

  matched_jet_pairs = []

  num_offlineEvents = 0

  # main loop that goes over the file
  for event_num in range(total_num_events):
    if event_num % 100 == 0 and event_num != 0:
      print "doing event_num=%d for (%d, %d, %d)" % (event_num, offline_jetpT_threshold, gTower_jetET_threshold, seed_ETthresh)
    # pull in data row by row
    data = rnp.root2rec(filename, treename='%s/%s' % (directory,tree), branches=offline_column_names + gTower_column_names, start=(event_num), stop=(event_num+1))
    oEvent = OfflineJets.Event(event=[data[col][0] for col in offline_column_names])

    # if there are no offline jets, we skip it
    if len(oEvent.jets) == 0 or oEvent.jets[0].pT < offline_jetpT_threshold:
      continue
    num_offlineEvents += 1

    '''can use seed_filter on an event by event basis'''
    # max number of seeds based on number of offline jets
    #seed_filter = gTowers.SeedFilter(numSeeds = len(oEvent.jets))
    tEvent = gTowers.TowerEvent(event=[data[col][0] for col in gTower_column_names], seed_filter = seed_filter)
    tEvent.get_event()

    #paired_jets = match_jets(oJets=oEvent.jets, tJets=tEvent.filter_towers())
    paired_jets = match_jets(oJets=oEvent.jets, tJets=tEvent.event.jets)
    matched_jet_pairs.append(np.array([[oJet.pT, tJet.E/np.cosh(tJet.eta)] for oJet,tJet in paired_jets if oJet.pT > offline_jetpT_threshold and tJet.E > 0. and tJet.E > 0.]))


  '''at this point, we've processed all the data and we just need to make plots'''

  filename_ending = 'offline%d_gTower%d_seed%d_unweighted' % (offline_jetpT_threshold, gTower_jetET_threshold, seed_filter.ETthresh)

  matched_jet_pairs = np.array(matched_jet_pairs)
  all_jet_pairs = np.array([l for item in matched_jet_pairs for l in item])

  leading_offline_jet_pairs = np.array([l for item in matched_jet_pairs for l in item if l[1] == np.amax(item[:,1])])

  xlim = (1e0,5e3)
  ylim = (1e1,1e4)

  #make figures
  '''All Jet Pairs'''
  pl_aJet = pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('trigger $E_T^{\mathrm{jet}}$ [GeV]')
  pl.title('$p_T^{\mathrm{offline jet}}$ > %d GeV, %d events, $E_T^{\mathrm{tower jet}}$ > %d GeV, $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.scatter(all_jet_pairs[:,0], all_jet_pairs[:,1])
  pl.xscale('log')
  pl.yscale('log')
  pl.grid(True, which='both')
  pl.xlim(xlim)
  pl.ylim(ylim)
  pickle.dump(pl_aJet, file('events_all_jet_pairs_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('events_all_jet_pairs_%s.png' % filename_ending)
  pl.close()

  '''Leading Offline Jet Pairs'''
  pl_lJet = pl.figure()
  pl.xlabel('leading offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('trigger $E_T^{\mathrm{jet}}$ [GeV]')
  pl.title('$p_T^{\mathrm{offline jet}}$ > %d GeV, %d events, $E_T^{\mathrm{tower jet}}$ > %d GeV, $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.scatter(leading_offline_jet_pairs[:,0], leading_offline_jet_pairs[:,1])
  pl.xscale('log')
  pl.yscale('log')
  pl.grid(True, which='both')
  pl.xlim(xlim)
  pl.ylim(ylim)
  pickle.dump(pl_lJet, file('events_leading_offline_jet_pairs_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('events_leading_offline_jet_pairs_%s.png' % filename_ending)
  pl.close()


class Copier(object):
  def __init__(self, offline_jetpT_threshold, gTower_jetET_threshold):
    self.offline_jetpT_threshold = offline_jetpT_threshold
    self.gTower_jetET_threshold = gTower_jetET_threshold
  def __call__(self, seed_ETthresh):
    run_code(self.offline_jetpT_threshold, self.gTower_jetET_threshold, seed_ETthresh)

import multiprocessing
p = multiprocessing.Pool(processes=6)
p.map(Copier(0., 0.), [10., 15., 20., 25., 30., 35.])
