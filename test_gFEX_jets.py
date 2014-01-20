from gFEX_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim.root'
tEvents = TowerEvents(filename=filename)
tEvents.load()

print 'Loaded and built all objects for gTowers data'
num = 0
grids = []
for tower_event in tEvents:
  grid = Grid(pixel_resolution=0.05)
  num += 1
  grid.add_tower_event(tower_event)
  grids.append(grid)
