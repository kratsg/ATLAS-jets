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

offline_column_names = ['jet_AntiKt10LCTopo_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi']]
gTower_column_names = ['gTower%s' % col for col in ['E', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]

seed_filter = gTowers.SeedFilter(20.0)

domain = np.array([[-3.2,3.2],[-4.9,4.9]])

def run_code(event_num):
  data = rnp.root2rec(filename, treename='%s/%s' % (directory,tree), branches=offline_column_names + gTower_column_names, start=(event_num), stop=(event_num+1))
  oEvent = OfflineJets.Event(event=[data[col][0] for col in offline_column_names])
  tEvent = gTowers.TowerEvent(event=[data[col][0] for col in gTower_column_names], seed_filter = seed_filter)
  tEvent.get_event()

  grid_towers = gTowers.Grid(cell_resolution=0.02, domain=domain)
  grid_towers.add_tower_event(tEvent)
  grid_towers.save(title='Event %d, gTowers, cell resolution=0.02' % event_num, filename='event_%d_towers.png' % event_num, colzLabel = '$E_T^{\mathrm{tower}}$')

  grid_offline = gTowers.Grid(cell_resolution=0.02, recon_algo = 'gaussian', domain=domain)
  grid_offline.add_event(oEvent)
  grid_offline.save(title='Event %d, offline jets, cell resolution=0.02' % event_num, filename='event_%d_offline_jets.png' % event_num, colzLabel = '$p_T^{\mathrm{jet}}$')

  grid_trigger = gTowers.Grid(cell_resolution=0.02, recon_algo = 'gaussian', domain=domain)
  grid_trigger.add_event(tEvent.get_event())
  grid_trigger.save(title='Event %d, trigger jets, cell resolution=0.02' % event_num, filename='event_%d_trigger_jets.png' % event_num, colzLabel = '$E_T^{\mathrm{jet}}$')

import multiprocessing
p = multiprocessing.Pool(processes=6)
p.map(run_code, range(total_num_events))
