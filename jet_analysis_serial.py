from atlas_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim_TTbar_14TeV_MU80.root'
directory = 'TTbar_14TeV_MU80/'
tree = 'mytree'

oEvents = OfflineJets.Events(rootfile)

grid = OfflineJets.Grid(pixel_resolution=0.2)
for oEvent in oEvents:
  grid.add_event(oEvent)

# at this point, all of the events have been processed (including jets)
analysis = Analysis(offline_events = oEvents)
analysis.Efficiency()
