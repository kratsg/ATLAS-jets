from gFEX_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim.root'
tEvents = TowerEvents(filename=filename)
tEvents.load()

print 'Loaded and built all objects for gTowers data'
