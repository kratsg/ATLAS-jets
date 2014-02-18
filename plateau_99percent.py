from atlas_jets import *
import ROOT
import numpy as np
import matplotlib.pyplot as pl
import root_numpy as rnp
import cPickle as pickle
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

bins_seedcut_99percent = []
bins_seedcut_errors = []

for seed_ETthresh in [5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0]:
  #set seed cuts
  seed_filter = gTowers.SeedFilter(ETthresh = seed_ETthresh, numSeeds = 1.0e5)

  #column names to pull from the file, must be in this order to sync with the predefined classes in atlas_jets package
  offline_column_names = ['jet_AntiKt10LCTopo_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi']]
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

    #tEvent.get_event()
    paired_jets = match_jets(oJets=oEvent.jets, tJets=tEvent.filter_towers())
    #paired_jets = match_jets(oJets=oEvent.jets, tJets=tEvent.event.jets)
    paired_data = np.array([[oJet.pT, tJet.E/np.cosh(tJet.eta)] for oJet,tJet in paired_jets if oJet.pT > offline_jetpT_threshold])
    # build up the turn on curve histograms
    hist_efficiency_den += np.histogram(paired_data[:,0], bins=bins_efficiency)[0]
    hist_efficiency_num += np.histogram(paired_data[np.where(paired_data[:,1] > gTower_jetET_threshold),0], bins=bins_efficiency)[0]


  '''at this point, we've processed all the data and we just need to make plots'''

  # first get the widths of the bins when we make the plots
  width_efficiency        = np.array([x - bins_efficiency[i-1] for i,x in enumerate(bins_efficiency)][1:])

  filename_ending = 'offline%d_gTower%d_seed%d_unweighted' % (offline_jetpT_threshold, gTower_jetET_threshold, seed_filter.ETthresh)

  xlim_efficiency = (0.0,1.0)
  ylim_efficiency = (0.0,1.0)

  '''Turn on curves'''
  pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Turn-On Curve Denominator')
  pl.title('offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events, gTower $E_T^{\mathrm{jet}}$ > %d GeV, gTower $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
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
  pickle.dump(pl_turnon_den, file('events_turnon_denominator_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('events_turnon_denominator_%s.png' % filename_ending)
  pl.close()

  pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Turn-On Curve Numerator')
  pl.title('offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events, gTower $E_T^{\mathrm{jet}}$ > %d GeV, gTower $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
  pl.bar(bins_efficiency[:-1], hist_efficiency_num, width=width_efficiency)
  pl.xlim(xlim_efficiency)
  pl.ylim(ylim_efficiency)
  pl_turnon_num = {'bins': bins_efficiency,\
                   'values': hist_efficiency_num,\
                   'width': width_efficiency}
  pickle.dump(pl_turnon_num, file('events_turnon_numerator_%s.pkl' % filename_ending, 'w+') )
  pl.savefig('events_turnon_numerator_%s.png' % filename_ending)
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
        errors.append( (np.abs( (1.-2.*w + w**2.)/den))**0.5 )
    return errors

  #binomial errors s^2 = n * p * q
  errors_efficiency_differential = binomial_errors(hist_efficiency_curve_differential, hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
  errors_efficiency_integral     = binomial_errors(hist_efficiency_curve_integral, np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])

  pl.figure()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Trigger Efficiency - Differential')
  pl.title('offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events, gTower $E_T^{\mathrm{jet}}$ > %d GeV, gTower $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
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
  pl.ylabel('Trigger Efficiency - Integral')
  pl.title('offline $p_T^{\mathrm{jet}}$ > %d GeV, %d events, gTower $E_T^{\mathrm{jet}}$ > %d GeV, gTower $E_T^{\mathrm{seed}}$ > %d GeV' % (offline_jetpT_threshold, num_offlineEvents, gTower_jetET_threshold, seed_filter.ETthresh))
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
  pickle.dump(pl_eff_int, file('events_turnon_curve_integral_%s.pkl' % filename_ending, 'w+'))
  pl.savefig('events_turnon_curve_integral_%s.png' % filename_ending)
  pl.close()

  xpoints_efficiency = xpoints_efficiency[nonzero_bins]
  index_99percent = np.where(hist_efficiency_curve_differential >= 0.99)[0][0]
  bins_around_99percent_eff = xpoints_efficiency[index_99percent-1:index_99percent+2]
  bins_seedcut_99percent.append([ seed_ETthresh,bins_around_99percent_eff[1] ])
  bins_seedcut_errors.append(bins_around_99percent_eff[2] - bins_around_99percent_eff[0])

bins_seedcut_99percent = np.array(bins_seedcut_99percent)
bins_seedcut_errors = np.array(bins_seedcut_errors)
pl.figure()
pl.xlabel('gTower $E_T^{\mathrm{seed}}$ threshold [GeV]')
pl.ylabel('offline $p_T^{\mathrm{jet}}$ [GeV]')
pl.title('99% Plateau with offline $p_T^{\mathrm{jet}}$ > 0 GeV, gTower $E_T^{\mathrm{jet}}$ > 0 GeV')
pl.errorbar(bins_seedcut_99percent[:,0], bins_seedcut_99percent[:,1], xerr=1.0, yerr=bins_seedcut_errors)
pl_plateau = {'xdata': bins_seedcut_99percent[:,0],\
              'ydata': bins_seedcut_99percent[:,1],\
              'xerr' : 1.0,\
              'yerr' : bins_seedcut_errors}
pl.grid(True)
pickle.dump(pl_plateau, file('events_plateau_99percent.pkl', 'w+'))
pl.savefig('events_plateau_99percent.png')

