from atlas_jets import *
import root_numpy as rnp
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim_TTbar_14TeV_MU80.root'
directory = 'TTbar_14TeV_MU80'
tree = 'mytree'

#set total number of events
total_num_events = 1000
#set the seed filter
seed_filter = gTowers.SeedFilter(ETthresh = 0.0, numSeeds = 10)


offline_column_names = ['jet_AntiKt4LCTopo_%s' % col for col in ['E', 'pt', 'm', 'eta', 'phi']]
gTower_column_names = ['gTower%s' % col for col in ['E', 'NCells', 'EtaMin', 'EtaMax', 'PhiMin', 'PhiMax']]

for event_num in range(total_num_events):
  # pull in data row by row
  data = rnp.root2rec(filename, treename='%s/%s' % (directory,tree), branches=offline_column_names + gTower_column_names, start=(event_num), stop=(event_num+1))
  oEvent = OfflineJets.Event(event=[data[col][0] for col in offline_column_names])
  tEvent = gTowers.TowerEvent(event=[data[col][0] for col in gTower_column_names], seed_filter = seed_filter)
  tEvent.get_event()
  print oEvent
  print tEvent
  print tEvent.event
