from atlas_jets import *
import ROOT
import numpy as np
import matplotlib.pyplot as pl
import root_numpy as rnp
import pickle
filename = '~/Desktop/PileupSkim_TTbar_14TeV_MU80_10000.root' 
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
  return np.array(matched_jets)

def run_code(offline_jetpT_threshold = 0., gTower_jetET_threshold = 0., seed_ETthresh = 0.):
  #set seed cuts
  seed_filter = gTowers.SeedFilter(ETthresh = seed_ETthresh, numSeeds = 1.0e5)

  leading_jet_data = []

  #column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
  offline_column_names = ['jet_AntiKt10LCTopo_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi']]
  gTower_column_names = ['gTower%s' % col for col in ['E', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]

  #bins for all histograms
  num_offlineEvents = 0

  # main loop that goes over the file
  for event_num in range(total_num_events):
    if event_num % 100 == 0 and event_num != 0:
      print "doing event_num=%d for (%d, %d, %d)" % (event_num, offline_jetpT_threshold, gTower_jetET_threshold, seed_ETthresh)
      break
    # pull in data row by row
    data = rnp.root2rec(filename, treename='%s/%s' % (directory,tree), branches=offline_column_names + gTower_column_names, start=(event_num), stop=(event_num+1))
    oEvent = OfflineJets.Event(event=[data[col][0] for col in offline_column_names])

    # if there are no offline jets, we skip it
    if len(oEvent.jets) == 0 or oEvent.jets[0].pT < offline_jetpT_threshold or np.abs(oEvent.jets[0].rapidity) > 1.2:
      continue
    num_offlineEvents += 1

    '''can use seed_filter on an event by event basis'''
    # max number of seeds based on number of offline jets
    #seed_filter = gTowers.SeedFilter(numSeeds = len(oEvent.jets))
    tEvent = gTowers.TowerEvent(event=[data[col][0] for col in gTower_column_names], seed_filter = seed_filter)
    # build up the first two histograms using just the gTower data
    # note, we have np.histogram(...)[0] since we only need the hist data
    denominator.append(oEvent.jets[0].pT)
    seeds = [seed for seed in tEvent.filter_towers() if np.abs(seed.eta) < 1.6]
    if len(seeds) > 0:
      numerator.append(oEvent.jets[0].pT)

  '''at this point, we've processed all the data and we just need to make plots'''

  # first get the widths of the bins when we make the plots
  bins_efficiency        = np.arange(0,1240, 20).astype(float)
  width_efficiency        = np.array([x - bins_efficiency[i-1] for i,x in enumerate(bins_efficiency)][1:])

  hist_efficiency_den, _ = np.histogram(np.array(denominator), bins=bins_efficiency)
  hist_efficiency_num, _ = np.histogram(np.array(numerator), bins=bins_efficiency)

  filename_ending = 'offline%d_gTower%d_seed%d_unweighted' % (offline_jetpT_threshold, gTower_jetET_threshold, seed_filter.ETthresh)

  '''Turn on curves'''
  pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Turn-On Curve Denominator')
  pl.title('$p_T^{\mathrm{offline jet}}$ > %d GeV, %d events, $E_T^{\mathrm{tower jet}}$ > %d GeV, $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.bar(bins_efficiency[:-1], hist_efficiency_den, width=width_efficiency)
  xlim_efficiency = pl.xlim()
  xlim_efficiency = (0.0, xlim_efficiency[1])
  pl.xlim(xlim_efficiency)
  ylim_efficiency = pl.ylim()
  ylim_efficiency = (0.0, ylim_efficiency[1])
  pl.ylim(ylim_efficiency)
  pl_turnon_den = {'bins': bins_efficiency,\
                   'values': hist_efficiency_den,\
                   'width': width_efficiency}
  pickle.dump(pl_turnon_den, file('eventwise_seed_turnon_denominator_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('eventwise_seed_turnon_denominator_%s.png' % filename_ending)
  pl.close()

  pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Turn-On Curve Numerator')
  pl.title('$p_T^{\mathrm{offline jet}}$ > %d GeV, %d events, $E_T^{\mathrm{tower jet}}$ > %d GeV, $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.bar(bins_efficiency[:-1], hist_efficiency_num, width=width_efficiency)
  pl.xlim(xlim_efficiency)
  pl.ylim(ylim_efficiency)
  pl_turnon_num = {'bins': bins_efficiency,\
                   'values': hist_efficiency_num,\
                   'width': width_efficiency}
  pickle.dump(pl_turnon_num, file('eventwise_seed_turnon_numerator_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('eventwise_seed_turnon_numerator_%s.png' % filename_ending)
  pl.close()

  nonzero_bins = np.where(hist_efficiency_den != 0)
  #compute integral and differential curves
  hist_efficiency_curve_differential = np.true_divide(hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
  hist_efficiency_curve_integral = np.true_divide(np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])
  #get halfway in between really
  xpoints_efficiency = bins_efficiency[:-1] + width_efficiency/2.

  def binomial_errors(hist_ratio, hist_one, hist_two):
    errors = []
    for w, num, den in zip(hist_ratio, hist_one, hist_two):
      # root.cern.ch/root/html/src/TH1.cxx.html#l5.yxD
      # formula cited (for histograms [num, den] with no errors) is:
      #     w = num/den
      #     if w = 1:
      #             sigma = 0
      #     else:
      #             sigma = abs( (1 - 2*w + w**2) / den**2 )
      if w == 1.0:
        errors.append(0.0)
      else:
        errors.append( (np.abs( (1.-2.*w + w**2.)/den**2.))**0.5 )
    return errors

  #binomial errors s^2 = n * p * q
  errors_efficiency_differential = binomial_errors(hist_efficiency_curve_differential, hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
  errors_efficiency_integral     = binomial_errors(hist_efficiency_curve_integral, np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])

  pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Seed Trigger Efficiency - Differential')
  pl.title('$p_T^{\mathrm{offline jet}}$ > %d GeV, %d events, $E_T^{\mathrm{tower jet}}$ > %d GeV, $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.errorbar(xpoints_efficiency[nonzero_bins], hist_efficiency_curve_differential, yerr=errors_efficiency_differential, ecolor='black')
  pl.xlim(xlim_efficiency)
  pl.ylim((0.0,1.2))
  pl.grid(True)
  pl_eff_diff = {'xdata': xpoints_efficiency,\
                 'ydata': hist_efficiency_curve_differential,\
                 'xerr' : 1.0,\
                 'yerr' : errors_efficiency_differential,\
                 'num'  : hist_efficiency_num,\
                 'den'  : hist_efficiency_den,\
                 'bins' : bins_efficiency,\
                 'nonzero_bins': nonzero_bins}
  pickle.dump(pl_eff_diff, file('events_turnon_curve_differential_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('events_turnon_curve_differential_%s.png' % filename_ending)
  pl.close()

  pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Seed Trigger Efficiency - Integral')
  pl.title('$p_T^{\mathrm{offline jet}}$ > %d GeV, %d events, $E_T^{\mathrm{tower jet}}$ > %d GeV, $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.errorbar(xpoints_efficiency[nonzero_bins], hist_efficiency_curve_integral, yerr=errors_efficiency_integral, ecolor='black')
  pl.xlim(xlim_efficiency)
  pl.ylim((0.0,1.2))
  pl.grid(True)
  pl_eff_int = {'xdata': xpoints_efficiency,\
                'ydata': hist_efficiency_curve_integral,\
                'xerr' : 1.0,\
                'yerr' : errors_efficiency_integral,\
                'num'  : hist_efficiency_num,\
                'den'  : hist_efficiency_den,\
                'bins' : bins_efficiency,\
                'nonzero_bins': nonzero_bins}
  pickle.dump(pl_eff_int, file('events_turnon_curve_integral_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('events_turnon_curve_integral_%s.png' % filename_ending)
  pl.close()

class Copier(object):
  def __init__(self):
    pass
  def __call__(self, args):
    offline_jetpT_threshold, gTower_jetET_threshold, seed_ETthresh = args
    run_code(offline_jetpT_threshold, gTower_jetET_threshold, seed_ETthresh)

offline_jetpT_threshold = [0.]
gTower_jetET_threshold = [0.]
seed_ETthreshold = [0.]
#seed_ETthreshold = [35., 30., 25., 20.]
jobs = [(a,b,c) for a in offline_jetpT_threshold for b in gTower_jetET_threshold for c in seed_ETthreshold]
print jobs

import multiprocessing
p = multiprocessing.Pool(processes=6)
p.map(Copier(), jobs)
