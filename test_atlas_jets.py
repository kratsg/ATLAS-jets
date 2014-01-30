from atlas_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim_TTbar_14TeV_MU80.root'
directory = 'TTbar_14TeV_MU80/'
tree = 'mytree'

rootfile = RootFile(filename=filename, directory=directory, tree=tree)
oEvents = OfflineJets.Events(rootfile)
tEvents = gTowers.TowerEvents(rootfile)
analysis = Analysis(oEvents, tEvents)

for pT_thresh in [0.,50.,100.,150.,200.,250.,300.,350.,400.]:
  analysis.TowerMultiplicity(pT_thresh = pT_thresh)
  analysis.TowerHistogram(pT_thresh = pT_thresh)


raise SystemExit
#this goes through and plots each event using the towers
for tEvent in tEvents:
  num += 1
  grid_before = gTowers.Grid(pixel_resolution=0.01, domain=domain)
  grid_before.add_tower_event(tEvent)
  grid_before.save(title='Event %d, cell resolution=0.01' % num, filename='output_tower_%d.png' % num, colzLabel = '$E_T^{\mathrm{tower}}$')

  grid_after = gTowers.Grid(pixel_resolution=0.01, recon_algo = 'gaussian', domain=domain)
  grid_after.add_event(tEvent.get_event())
  grid_after.save(title='Event %d, cell resolution=0.01' % num, filename='output_event_%d.png' % num, colzLabel = '$p_T^{\mathrm{jet}}$')

