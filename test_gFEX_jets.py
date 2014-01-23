from gFEX_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim.root'
tEvents = TowerEvents(filename=filename)
tEvents.load()

print 'Loaded and built all objects for gTowers data'
num = 0
grids = []

#this goes through and plots each event using the towers
for tower_event in tEvents:
  break #remoe if you want to go through loop
  grid = Grid(pixel_resolution=0.01)
  num += 1
  grid.add_tower_event(tower_event)
  grids.append(grid)
  grid.save(title='Event %d, cell resolution=0.01' % num, filename='output_event_%d.png' % num)

seed_filter = SeedFilter(ETthresh = 20.0)#20.0 GeV
for tEvent in tEvents:
  tEvent.set_seed_filter(seed_filter)

