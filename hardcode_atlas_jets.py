from atlas_jets import *
import ROOT
import numpy as np
import matplotlib.pyplot as pl
import root_numpy as rnp
filename = '/Users/kratsg/Desktop/PileupSkim_TTbar_14TeV_MU80_10000.root'
directory = 'TTbar_14TeV_MU80'
tree = 'mytree'

rootfile = ROOT.TFile(filename)
#set total number of events
total_num_events = int(rootfile.Get('%s/%s' % (directory, tree)).GetEntries())

#set jet cuts
offline_jetpT_threshold = 200. #[GeV]
gTower_jetET_threshold  = 150.

#set seed cuts
seed_filter = gTowers.SeedFilter(ETthresh = 20.0, numSeeds = 1.0e5)

#column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
offline_column_names = ['jet_AntiKt4LCTopo_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi']]
gTower_column_names = ['gTower%s' % col for col in ['E', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]

#bins for all histograms
bins_towerMultiplicity = np.arange(0, 1000, 5).astype(float)
bins_towerHistogram    = np.array([0,50,100,150,200,250,300,350,400,500,750,1000,4000]).astype(float)
bins_efficiency        = np.arange(0,1240, 20).astype(float)

hist_towerMultiplicity = np.zeros(len(bins_towerMultiplicity)-1).astype(float)
hist_towerHistogram    = np.zeros(len(bins_towerHistogram)-1).astype(float)
hist_efficiency_num    = np.zeros(len(bins_efficiency)-1).astype(float)
hist_efficiency_den    = np.zeros(len(bins_efficiency)-1).astype(float)

num_offlineEvents = 0

# define helper functions - also a source of parallelization!
def compute_jetDistance(jet1, jet2):
  return ((jet1.eta - jet2.eta)**2. + (jet1.phi - jet2.phi)**2.)**0.5

def match_jets(oJets=[], tJets=[]):
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

# main loop that goes over the file
for event_num in range(total_num_events):
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
  # build up the first two histograms using just the gTower data
  # note, we have np.histogram(...)[0] since we only need the hist data
  tower_ETs = [tower.E/np.cosh(tower.eta) for tower in tEvent.towers]
  hist_towerMultiplicity += np.cumsum(np.histogram(tower_ETs, bins=bins_towerMultiplicity)[0][::-1])[::-1] #this makes a reverse cumulative sum
  hist_towerHistogram += np.histogram(tower_ETs, bins=bins_towerHistogram)[0]

  tEvent.get_event()
  paired_jets = match_jets(oJets=oEvent.jets, tJets=tEvent.event.jets)
  paired_data = np.array([[oJet.pT, tJet.E/np.cosh(tJet.eta)] for oJet,tJet in paired_jets if oJet.pT > offline_jetpT_threshold])
  # build up the turn on curve histograms
  hist_efficiency_den += np.histogram(paired_data[:,0], bins=bins_efficiency)[0]
  hist_efficiency_num += np.histogram(paired_data[np.where(paired_data[:,1] > gTower_jetET_threshold),0], bins=bins_efficiency)[0]


'''at this point, we've processed all the data and we just need to make plots'''

# first get the widths of the bins when we make the plots
width_towerMultiplicity = np.array([x - bins_towerMultiplicity[i-1] for i,x in enumerate(bins_towerMultiplicity)][1:])
width_towerHistogram    = np.array([x - bins_towerHistogram[i-1] for i,x in enumerate(bins_towerHistogram)][1:])
width_efficiency        = np.array([x - bins_efficiency[i-1] for i,x in enumerate(bins_efficiency)][1:])

# rescale tower data to define it per event
hist_towerMultiplicity = 1.0*hist_towerMultiplicity/num_offlineEvents
hist_towerHistogram    = 1.0*hist_towerHistogram/num_offlineEvents

#histogram y-range
hist_ylim = (10.**-3., 10.**4.)

#make figures
'''Tower Multiplicity'''
pl.figure()
pl.xlabel('$E_T^{\mathrm{threshold}}$ [GeV]')
pl.ylabel('Number of gTowers per event')
pl.title('Number of gTowers above $E_T^{\mathrm{threshold}}$ for offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events' % (offline_jetpT_threshold, num_offlineEvents))
pl.bar(bins_towerMultiplicity[:-1], hist_towerMultiplicity, width=width_towerMultiplicity, log=True)
pl.ylim(hist_ylim)
pl.savefig('events_threshold_histogram_multiplicity_offline%d_gTower%d.png' % (offline_jetpT_threshold, gTower_jetET_threshold))
pl.close()

'''Tower Histogram per Event'''
pl.figure()
pl.xlabel('$p_T^{\mathrm{jet}}$ [GeV]')
pl.ylabel('Number of gTowers per event')
pl.title('Histogram of gTowers for offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events' % (offline_jetpT_threshold, num_offlineEvents))
pl.bar(bins_towerHistogram[:-1], hist_towerHistogram, width=width_towerHistogram, log=True)
pl.ylim(hist_ylim)
pl.savefig('events_threshold_histogram_towers_offline%d_gTower%d.png' % (offline_jetpT_threshold, gTower_jetET_threshold))
pl.close()

xlim_efficiency = (0.0,1.0)
ylim_efficiency = (0.0,1.0)

'''Turn on curves'''
pl.figure()
pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
pl.ylabel('Efficiency')
pl.title('Turn-On Denominator for offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events' % (offline_jetpT_threshold, num_offlineEvents))
pl.bar(bins_efficiency[:-1], hist_efficiency_den, width=width_efficiency)
xlim_efficiency = pl.xlim()
xlim_efficiency = (0.0, xlim_efficiency[1])
pl.xlim(xlim_efficiency)
ylim_efficiency = pl.ylim()
ylim_efficiency = (0.0, ylim_efficiency[1])
pl.ylim(ylim_efficiency)
pl.savefig('events_turnon_denominator_offline%d_gTower%d.png' % (offline_jetpT_threshold, gTower_jetET_threshold))
pl.close()

pl.figure()
pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
pl.ylabel('Efficiency')
pl.title('Turn-On Numerator for offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events' % (offline_jetpT_threshold, num_offlineEvents))
pl.bar(bins_efficiency[:-1], hist_efficiency_num, width=width_efficiency)
pl.xlim(xlim_efficiency)
pl.ylim(ylim_efficiency)
pl.savefig('events_turnon_numerator_offline%d_gTower%d.png' % (offline_jetpT_threshold, gTower_jetET_threshold))
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
pl.ylabel('Trigger Efficiency')
pl.title('Differential Turn-On Curve for offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events' % (offline_jetpT_threshold, num_offlineEvents))
pl.errorbar(xpoints_efficiency[nonzero_bins], hist_efficiency_curve_differential, yerr=errors_efficiency_differential, ecolor='black')
pl.xlim(xlim_efficiency)
pl.ylim((0.0,1.2))
pl.grid(True)
pl.savefig('events_turnon_curve_differential_offline%d_gTower%d.png' % (offline_jetpT_threshold, gTower_jetET_threshold))
pl.close()

pl.figure()
pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
pl.ylabel('Trigger Efficiency')
pl.title('Integral Turn-On Curve for offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events' % (offline_jetpT_threshold, num_offlineEvents))
pl.errorbar(xpoints_efficiency[nonzero_bins], hist_efficiency_curve_integral, yerr=errors_efficiency_integral, ecolor='black')
pl.xlim(xlim_efficiency)
pl.ylim((0.0,1.2))
pl.grid(True)
pl.savefig('events_turnon_curve_integral_offline%d_gTower%d.png' % (offline_jetpT_threshold, gTower_jetET_threshold))
pl.close()

