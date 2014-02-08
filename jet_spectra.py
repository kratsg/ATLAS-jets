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

def run_code(offline_jetpT_threshold = 0., gTower_jetET_threshold = 0., seed_ETthresh = 0.):
  #set seed cuts
  seed_filter = gTowers.SeedFilter(ETthresh = seed_ETthresh, numSeeds = 1)

  leading_trigger_jets = []

  #column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
  offline_column_names = ['jet_AntiKt4LCTopo_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi']]
  gTower_column_names = ['gTower%s' % col for col in ['E', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]

  #bins for all histograms
  num_offlineEvents = 0

  # main loop that goes over the file
  for event_num in range(total_num_events):
    if event_num % 100 == 0:
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
    leading_trigger_jets = leading_trigger_jets + tEvent.event.jets


  '''at this point, we've processed all the data and we just need to make plots'''

  bins_leading_trigger_jets  = np.arange(0.,4000.,5.)
  hist_leading_trigger_jets  = np.histogram([jet.E/np.cosh(jet.eta) for jet in leading_trigger_jets], bins=bins_leading_trigger_jets)[0]
  # first get the widths of the bins when we make the plots
  width_leading_trigger_jets = np.array([x - bins_leading_trigger_jets[i-1] for i,x in enumerate(bins_leading_trigger_jets)][1:])

  filename_ending = 'offline%d_gTower%d_seed%d_unweighted' % (offline_jetpT_threshold, gTower_jetET_threshold, seed_filter.ETthresh)

  #make figures
  '''Leading Trigger Jets Histogram'''
  pl_lJet = pl.figure()
  pl.xlabel('$E_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Number of leading trigger jets')
  pl.title('$p_T^{\mathrm{offline jet}}$ > %d GeV, %d events, $E_T^{\mathrm{tower jet}}$ > %d GeV, $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.bar(bins_leading_trigger_jets[:-1], hist_leading_trigger_jets, width=width_leading_trigger_jets)
  pickle.dump(pl_lJet, file('events_histogram_leading_trigger_jets_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('events_histogram_leading_trigger_jets_%s.png' % filename_ending)
  pl.close()


class Copier(object):
  def __init__(self, offline_jetpT_threshold, gTower_jetET_threshold):
    self.offline_jetpT_threshold = offline_jetpT_threshold
    self.gTower_jetET_threshold = gTower_jetET_threshold
  def __call__(self, seed_ETthresh):
    run_code(self.offline_jetpT_threshold, self.gTower_jetET_threshold, seed_ETthresh)

import multiprocessing
p = multiprocessing.Pool(processes=6)
p.map(Copier(0., 0.), [5., 10., 15., 20., 25., 30., 35., 40.])
