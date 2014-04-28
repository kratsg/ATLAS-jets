import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
import pickle
import re
from glob import glob
import os

p = re.compile('^(unmatched|matched)_jets_seed(\d+)_unweighted_page(?:\d+).pkl$')

bins_efficiency = np.arange(0.,2000.,10.)
width_efficiency = np.array([x - bins_efficiency[i-1] for i,x in enumerate(bins_efficiency)][1:])

files = glob("*matched_jets_seed*_unweighted_page*.pkl")

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

#this is a wrapper around file strings to ensure the directory exists
#       could use decorators...
def write_file(f):
  ensure_dir(f)
  return f

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
      errors.append( (np.abs( (1.-2.*w + w**2.)/den**2.))**0.5/2. )
  return errors

def distance(etas, phis):
  return ((etas[0] - etas[1])**2. + (phis[0] - phis[1])**2.)**0.5

def isolation_condition(jets):
  if np.abs(jets[0][0].eta) > 2.5:
    return False
  for oJet, tJet in jets:
    #skip leading jet
    if oJet == jets[0][0]:
      continue
    if distance([oJet.eta, jets[0][0].eta], [oJet.phi, jets[0][0].phi]) < 0.2:
      return False
  return True

figsize = (16, 12)
labelsize = 28
titlesize = 36

raw_data = []
for filename in files:
  matched, seed = p.match(filename).groups()
  seed = int(seed)
  step = '4'
  if step in ['3','4']:
    ylabel = "%s trigger $E_T^{\mathrm{jet}}$ [GeV]" % matched
  else:
    ylabel = "%s gTower $E_T^{\mathrm{tower}}$ [GeV]" % matched
  print "loading"
  raw_data += pickle.load(file(filename))
  print "loaded"

def nsj_filter(nsj, config):
  return eval(config)

pl.figure(figsize=figsize)
pl.hist([jets[0][0].nsj for jets in raw_data], bins=np.arange(0.,15.,1.))
pl.xlabel('Offline Jet - Number of Subjets', fontsize=labelsize)
pl.ylabel('count', fontsize=labelsize)
pl.title('Histogram of number of subjets for offline jets', fontsize=titlesize)
pl.savefig( write_file("TurnOnCurves/2dhistograms/step%s_nsj.png" % (step)) )
pl.close()

bins_split = (np.arange(0.,20.,1.),np.arange(0.,300000.,1000.))

pl.figure(figsize=figsize)
pl.hist2d([jets[0][0].nsj for jets in raw_data], [jets[0][0].split[0] for jets in raw_data], norm=LogNorm(), bins=bins_split )
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel('Offline Jet - Number of Subjets', fontsize=labelsize)
pl.ylabel('Offline Jet - Split12', fontsize=labelsize)
pl.title('Correlation of nsj and split12', fontsize=titlesize)
pl.xlim((0.,15.))
pl.ylim((0.,200000.))
pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_split12.png' % (step)) )
pl.close()

pl.figure(figsize=figsize)
pl.hist2d([jets[0][0].nsj for jets in raw_data], [jets[0][0].split[1] for jets in raw_data], norm=LogNorm(), bins=bins_split)
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel('Offline Jet - Number of Subjets', fontsize=labelsize)
pl.ylabel('Offline Jet - Split23', fontsize=labelsize)
pl.title('Correlation of nsj and split23', fontsize=titlesize)
pl.xlim((0.,15.))
pl.ylim((0.,200000.))
pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_split23.png' % (step)) )
pl.close()

pl.figure(figsize=figsize)
pl.hist2d([jets[0][0].nsj for jets in raw_data], [jets[0][0].split[2] for jets in raw_data], norm=LogNorm(), bins=bins_split)
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel('Offline Jet - Number of Subjets', fontsize=labelsize)
pl.ylabel('Offline Jet - Split34', fontsize=labelsize)
pl.title('Correlation of nsj and split34', fontsize=titlesize)
pl.xlim((0.,15.))
pl.ylim((0.,200000.))
pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_split34.png' % (step)) )
pl.close()

bins_tau = (np.arange(0.,20.,1.),np.arange(0.,2.,0.01))

pl.figure(figsize=figsize)
pl.hist2d([jets[0][0].nsj for jets in raw_data], [jets[0][0].tau[0] for jets in raw_data], norm=LogNorm(), bins=bins_tau)
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel('Offline Jet - Number of Subjets', fontsize=labelsize)
pl.ylabel('Offline Jet - $\\tau_1$', fontsize=labelsize)
pl.title('Correlation of nsj and $\\tau_1$', fontsize=titlesize)
pl.xlim((0.,15.))
pl.ylim((0.,1.))
pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_tau1.png' % (step)) )
pl.close()

pl.figure(figsize=figsize)
pl.hist2d([jets[0][0].nsj for jets in raw_data], [jets[0][0].tau[1] for jets in raw_data], norm=LogNorm(), bins=bins_tau)
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel('Offline Jet - Number of Subjets', fontsize=labelsize)
pl.ylabel('Offline Jet - $\\tau_2$', fontsize=labelsize)
pl.title('Correlation of nsj and $\\tau_2$', fontsize=titlesize)
pl.xlim((0.,15.))
pl.ylim((0.,1.))
pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_tau2.png' % (step)) )
pl.close()

pl.figure(figsize=figsize)
pl.hist2d([jets[0][0].nsj for jets in raw_data], [jets[0][0].tau[2] for jets in raw_data], norm=LogNorm(), bins=bins_tau)
cbar = pl.colorbar()
cbar.set_label('number density', fontsize=labelsize)
pl.xlabel('Offline Jet - Number of Subjets', fontsize=labelsize)
pl.ylabel('Offline Jet - $\\tau_3$', fontsize=labelsize)
pl.title('Correlation of nsj and $\\tau_3$', fontsize=titlesize)
pl.xlim((0.,15.))
pl.ylim((0.,1.))
pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_tau3.png' % (step)) )
pl.close()

# True == all nsj
for nsj_str, nsj_config in {'nsj(all)':'True', 'nsj(eq1)':'nsj == 1', 'nsj(eq2)':'nsj == 2', 'nsj(eq3)':'nsj == 3', 'nsj(eq4)':'nsj==4', 'nsj(geq2)':'nsj >= 2', 'nsj(geq3)':'nsj >= 3'}.iteritems():

  data = np.array([ [jets[0][0].Pt, jets[0][1].vector.Et()] for jets in raw_data if nsj_filter(jets[0][0].nsj, nsj_config) and isolation_condition(jets) ])

  pl.figure(figsize=figsize)
  pl.hist2d(data[:,0][np.where(data[:,1] > 0)], data[:,1][np.where(data[:,1] > 0)], bins=np.arange(0.,1500.,25.), norm=LogNorm())
  cbar = pl.colorbar()
  pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
  pl.ylabel(ylabel, fontsize=labelsize)
  pl.title('Correlation plot for Step %s, Seed Cut %d GeV' % (step, seed), fontsize=titlesize)
  cbar.set_label('number density', fontsize=labelsize)
  pl.xlim((0.0, 1500.0))
  pl.ylim((0.0, 1500.0))
  pl.savefig( write_file("TurnOnCurves/2dhistograms/step%s_seed%d_%s.png" %  (step, seed, nsj_str)) )
  pl.close()

  eta = np.array([[jets[0][0].eta, jets[0][1].eta] for jets in raw_data if jets[0][1].E != 0 and nsj_filter(jets[0][0].nsj, nsj_config) and isolation_condition(jets) ])
  phi = np.array([[jets[0][0].phi, jets[0][1].phi] for jets in raw_data if jets[0][1].E != 0 and nsj_filter(jets[0][0].nsj, nsj_config) and isolation_condition(jets) ])

  distances = np.array([distance(etas, phis) for etas, phis in zip(eta, phi)])

  pl.figure(figsize=figsize)
  pl.xlabel('Leading Offline $\eta^\mathrm{jet}$', fontsize=labelsize)
  pl.ylabel('Leading gTower $\eta^\mathrm{tower}$', fontsize=labelsize)
  pl.title('2D Correlation of $\eta$', fontsize=titlesize)
  pl.hist2d(eta[:,0],eta[:,1],bins=np.arange(-4.9,4.9,0.2), norm=LogNorm())
  pl.xlim((-4.9,4.9))
  pl.ylim((-4.9,4.9))
  cbar = pl.colorbar()
  cbar.set_label('number density', fontsize=labelsize)
  pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_seed%d_%s_eta.png' % (step, seed, nsj_str)) )
  pl.close()

  pl.figure(figsize=figsize)
  pl.xlabel('Leading Offline $\phi^\mathrm{jet}$', fontsize=labelsize)
  pl.ylabel('Leading gTower $\phi^\mathrm{tower}$', fontsize=labelsize)
  pl.title('2D Correlation of $\phi$', fontsize=titlesize)
  pl.hist2d(phi[:,0],phi[:,1],bins=np.arange(-4.9,4.9,0.2), norm=LogNorm())
  pl.xlim((-4.9,4.9))
  pl.ylim((-4.9,4.9))
  cbar = pl.colorbar()
  cbar.set_label('number density', fontsize=labelsize)
  pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_seed%d_%s_phi.png' % (step, seed, nsj_str)) )
  pl.close()

  pl.figure(figsize=figsize)
  pl.xlabel('Angular distance histogram', fontsize=labelsize)
  pl.ylabel('counts', fontsize=labelsize)
  pl.title('Histogram of angular distances', fontsize=titlesize)
  pl.hist(distances, bins=np.arange(0.0,1.0,0.05))
  pl.xlim((0.,1.))
  pl.savefig( write_file('TurnOnCurves/2dhistograms/step%s_seed%d_%s_angular_distance.png' % (step, seed, nsj_str)) )
  pl.close()

  # since these are for step 3 and step 4 (unmatched and matched)
  for trigger_jet_ETthresh in [0., 50., 100., 150., 200., 250., 300., 350., 400.]:
    hist_efficiency_den, _ = np.histogram(data[:,0], bins=bins_efficiency)
    hist_efficiency_num, _ = np.histogram(data[:,0][np.where(data[:,1] > trigger_jet_ETthresh)], bins=bins_efficiency)

    raw_energies = [[[oJet.Pt, tJet.vector.Et()] for oJet, tJet in jets] for jets in raw_data if nsj_filter(jets[0][0].nsj, nsj_config) and isolation_condition(jets) ]
    energies = np.array([[np.max([oJetPt for oJetPt, _ in event]),np.max([tJetEt for _, tJetEt in event])] for event in raw_energies])
        
    hist_J100_num, _ = np.histogram(energies[:,0][np.where(energies[:,1] > trigger_jet_ETthresh)], bins=bins_efficiency)

    pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Turn-On Curve Denominator', fontsize=labelsize)
    pl.title('Trigger Cut %d GeV, Seed Cut %d GeV, for Step %s' % (trigger_jet_ETthresh, seed, step), fontsize=titlesize)
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
    pickle.dump(pl_turnon_den, file( write_file('TurnOnCurves/Step%s/pickle/%s_turnon_denominator_step%s_seed%d_trigger%d_%s.pkl' % (step, matched, step, seed, trigger_jet_ETthresh, nsj_str)), 'w+') )
    pl.savefig( write_file('TurnOnCurves/Step%s/denominator/seed%d_trigger%d_%s.png' % (step, seed, trigger_jet_ETthresh, nsj_str)) )
    pl.close()

    pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Turn-On Curve Numerator', fontsize=labelsize)
    pl.title('Trigger Cut %d GeV, Seed Cut %d GeV for Step %s' % (trigger_jet_ETthresh, seed, step), fontsize=titlesize)
    pl.bar(bins_efficiency[:-1], hist_efficiency_num, width=width_efficiency)
    pl.xlim(xlim_efficiency)
    pl.ylim(ylim_efficiency)
    pl_turnon_num = {'bins': bins_efficiency,\
                     'values': hist_efficiency_num,\
                     'width': width_efficiency}
    pickle.dump(pl_turnon_num, file( write_file('TurnOnCurves/Step%s/pickle/%s_turnon_numerator_step%s_seed%d_trigger%d_%s.pkl' % (step, matched, step, seed, trigger_jet_ETthresh, nsj_str)), 'w+') )
    pl.savefig( write_file('TurnOnCurves/Step%s/numerator/seed%d_trigger%d_%s.png' % (step, seed, trigger_jet_ETthresh, nsj_str)) )
    pl.close()

    pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('J100 Turn-On Numerator', fontsize=labelsize)
    pl.title('Trigger Cut %d GeV, Seed Cut %d GeV for Step %s' % (trigger_jet_ETthresh, seed, step), fontsize=titlesize)
    pl.bar(bins_efficiency[:-1], hist_J100_num, width=width_efficiency)
    pl.xlim(xlim_efficiency)
    pl.ylim(ylim_efficiency)
    pl_J100_num = {'bins': bins_efficiency,\
                     'values': hist_J100_num,\
                     'width': width_efficiency}
    pickle.dump(pl_J100_num, file( write_file('TurnOnCurves/Step%s/pickle/%s_J100_numerator_step%s_seed%d_trigger%d_%s.pkl' % (step, matched, step, seed, trigger_jet_ETthresh, nsj_str)), 'w+') )
    pl.savefig( write_file('TurnOnCurves/Step%s/numerator/J100_seed%d_trigger%d_%s.png' % (step, seed, trigger_jet_ETthresh, nsj_str)) )
    pl.close()

    nonzero_bins = np.where(hist_efficiency_den != 0)
    #compute integral and differential curves
    hist_efficiency_curve_differential = np.true_divide(hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
    hist_efficiency_curve_integral = np.true_divide(np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])
    hist_J100_curve = np.true_divide(hist_J100_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
    #get halfway in between really
    xpoints_efficiency = bins_efficiency[:-1] + width_efficiency/2.

    #binomial errors s^2 = n * p * q
    errors_efficiency_differential = binomial_errors(hist_efficiency_curve_differential, hist_efficiency_num[nonzero_bins], hist_efficiency_den[nonzero_bins])
    errors_efficiency_integral     = binomial_errors(hist_efficiency_curve_integral, np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1], np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1])
    errors_J100_curve              = binomial_errors(hist_J100_curve, hist_J100_num[nonzero_bins], hist_efficiency_den[nonzero_bins])

    pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
    pl.title('Trigger Cut %d GeV, Seed Cut %d GeV, for Step %s' % (trigger_jet_ETthresh, seed, step), fontsize=titlesize)
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
    pickle.dump(pl_eff_diff, file( write_file('TurnOnCurves/Step%s/pickle/%s_turnon_differential_step%s_seed%d_trigger%d_%s.pkl' % (step, matched, step, seed, trigger_jet_ETthresh, nsj_str)), 'w+') )
    pl.savefig( write_file('TurnOnCurves/Step%s/differential/seed%d_trigger%d_%s.png' % (step, seed, trigger_jet_ETthresh, nsj_str)) )
    pl.close()

    pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Integral', fontsize=labelsize)
    pl.title('Trigger Cut %d GeV, Seed Cut %d GeV, for Step %s' % (trigger_jet_ETthresh, seed, step), fontsize=titlesize)
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
    pickle.dump(pl_eff_diff, file( write_file('TurnOnCurves/Step%s/pickle/%s_turnon_integral_step%s_seed%d_trigger%d_%s.pkl' % (step, matched, step, seed, trigger_jet_ETthresh, nsj_str)), 'w+') )
    pl.savefig( write_file('TurnOnCurves/Step%s/integral/seed%d_trigger%d_%s.png' % (step, seed, trigger_jet_ETthresh, nsj_str)) )
    pl.close()

    pl.figure(figsize=figsize)
    pl.xlabel('offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Event Efficiency - Integral', fontsize=labelsize)
    pl.title('Trigger Cut %d GeV, Seed Cut %d GeV, for Step %s' % (trigger_jet_ETthresh, seed, step), fontsize=titlesize)
    pl.errorbar(xpoints_efficiency[nonzero_bins], hist_J100_curve, yerr=errors_J100_curve, ecolor='black')
    pl.xlim(xlim_efficiency)
    pl.ylim((0.0,1.2))
    pl.grid(True)
    pl_J100_curve = {'xdata': xpoints_efficiency,\
                  'ydata': hist_J100_curve,\
                  'xerr' : 1.0,\
                  'yerr' : errors_J100_curve,\
                  'num'  : hist_J100_num,\
                  'den'  : hist_efficiency_den,\
                  'bins' : bins_efficiency,\
                  'nonzero_bins': nonzero_bins}
    pickle.dump(pl_J100_curve, file( write_file('TurnOnCurves/Step%s/pickle/%s_J100_curve_step%s_seed%d_trigger%d_%s.pkl' % (step, matched, step, seed, trigger_jet_ETthresh, nsj_str)), 'w+') )
    pl.savefig( write_file('TurnOnCurves/Step%s/J100/seed%d_trigger%d_%s.png' % (step, seed, trigger_jet_ETthresh, nsj_str)) )
    pl.close()
