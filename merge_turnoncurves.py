import numpy as np
import matplotlib.pyplot as pl
import pickle

steps = ['1','2','2prime','3','4']
seeds = [str(i) for i in range(0,60,5)]

def filename(step, seed, trigger = '0'):
  if step in ['1','3']:
    matched = "unmatched"
  else:
    matched = "matched"
  if step in ['1','2','2prime']:
    filename = "LeadingTurnOn/Step%s/pickle/%s_turnon_differential_step%s_seed%s.pkl" % (step, matched, step, seed)
  else:
    filename = "LeadingTurnOn/Step%s/pickle/%s_turnon_differential_step%s_seed%s_trigger%s.pkl" % (step, matched, step, seed, trigger)
  #print step, seed, "\n\t", filename
  return filename

figsize = (16,12)
ylim = (0.0,1.2)
xlim = (0.0,800.0)
labelsize = 28
titlesize = 48
legendprop = prop={'size':labelsize}
legendloc = 4

for step in steps:
  pl.figure(figsize=figsize)
  pl.xlabel('leading offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
  pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
  pl.title('Turn On Curves for Step %s' % step, fontsize=titlesize)
  pl.grid(True)
  for seed in seeds:
    if seed == '0' and step in ['3','4']:
      continue
    data = pickle.load(file( filename(step, seed) ))
    #xpoints_efficiency = data['xdata']
    #hist_efficiency_curve_differential = data['ydata']
    #errors_efficiency_differential = data['yerr']
    #nonzero_bins = data['nonzero_bins']
    pl.errorbar(data['xdata'][ data['nonzero_bins'] ], data['ydata'], yerr=data['yerr'], ecolor='black', label='Seed Cut %s GeV' % seed, linewidth=2)
  pl.ylim(ylim)
  pl.xlim(xlim)
  pl.legend(prop=legendprop, loc=legendloc)
  pl.savefig('MergedPlots/seeds/step%s.png' % step)
  pl.close()

markers = {'1': 'o', '2': 'v', '2prime': '^', '3': '<', '4': '>'}
for seed in seeds:
  pl.figure(figsize=figsize)
  pl.xlabel('leading offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
  pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
  pl.title('Turn On Curves for Seed Cut %s GeV' % seed, fontsize=titlesize)
  pl.grid(True)
  for step in steps:
    if seed == '0' and step in ['3','4']:
      continue
    data = pickle.load(file( filename(step, seed) ))
    pl.errorbar(data['xdata'][ data['nonzero_bins'] ], data['ydata'], yerr=data['yerr'], ecolor='black', label='Step %s' % step, linewidth=2, marker=markers[step], markersize=15)
  pl.ylim(ylim)
  pl.xlim(xlim)
  pl.legend(prop=legendprop, loc=legendloc)
  pl.savefig('MergedPlots/steps/seed%s.png' % seed)
  pl.close()

triggers = [str(i) for i in range(0,450,50)]

for seed in seeds[1:]:
  for step in ['3','4']:
    pl.figure(figsize=figsize)
    pl.xlabel('leading offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
    pl.title('Turn On Curves for Step %s, Seed Cut %s GeV' % (step, seed), fontsize=titlesize)
    pl.grid(True)
    for trigger in triggers:
      data = pickle.load(file( filename(step, seed, trigger) ))
      pl.errorbar(data['xdata'][ data['nonzero_bins'] ], data['ydata'], yerr=data['yerr'], ecolor='black', label='Trigger Cut %s GeV' % trigger, linewidth=2)
    pl.ylim(ylim)
    pl.xlim(xlim)
    pl.legend(prop=legendprop, loc=legendloc)
    pl.savefig('MergedPlots/triggers/step%s_seed%s.png' % (step, seed))
    pl.close()

for seed in seeds[1:]:
  for trigger in triggers:
    pl.figure(figsize=figsize)
    pl.xlabel('leading offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
    pl.title('Turn On Curves for Trigger Cut %s GeV, Seed Cut %s GeV' % (trigger, seed), fontsize=titlesize)
    pl.grid(True)
    for step in ['3','4']:
      data = pickle.load(file( filename(step, seed, trigger) ))
      pl.errorbar(data['xdata'][ data['nonzero_bins'] ], data['ydata'], yerr=data['yerr'], ecolor='black', label='Step %s' % step, linewidth=2, marker=markers[step], markersize=15)
    pl.ylim(ylim)
    pl.xlim(xlim)
    pl.legend(prop=legendprop, loc=legendloc)
    pl.savefig('MergedPlots/steps/trigger%s_seed%s.png' % (trigger, seed))
    pl.close()

for step in ['3','4']:
  for trigger in triggers:
    pl.figure(figsize=figsize)
    pl.xlabel('leading offline $p_T^{\mathrm{jet}}$ [GeV]', fontsize=labelsize)
    pl.ylabel('Trigger Efficiency - Differential', fontsize=labelsize)
    pl.title('Turn On Curves for Trigger %s GeV, Step %s' % (trigger, step), fontsize=titlesize)
    pl.grid(True)
    for seed in seeds[1:]:
      data = pickle.load(file( filename(step, seed, trigger) ))
      pl.errorbar(data['xdata'][ data['nonzero_bins'] ], data['ydata'], yerr=data['yerr'], ecolor='black', label='Seed Cut %s GeV' % seed, linewidth=2)
    pl.ylim(ylim)
    pl.xlim(xlim)
    pl.legend(prop=legendprop, loc=legendloc)
    pl.savefig('MergedPlots/seeds/step%s_trigger%s.png' % (step, trigger))
    pl.close()
