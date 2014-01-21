from gFEX_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim.root'
tEvents = TowerEvents(filename=filename)
tEvents.load()

print 'Loaded and built all objects for gTowers data'
num = 0
grids = []
for tower_event in tEvents:
  grid = Grid(pixel_resolution=0.01)
  num += 1
  grid.add_tower_event(tower_event)
  grids.append(grid)
  grid.save(title='Event %d, cell resolution=0.01' % num, filename='output_event_%d.png' % num)
