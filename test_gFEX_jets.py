from gFEX_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim.root'
tEvents = TowerEvents(filename=filename, seed_filter=SeedFilter(ETthresh=20.0))
tEvents.load()

print 'Loaded and built all objects for gTowers data'
num = 0

domain = np.array([[-3.2,3.2],[-4.9,4.9]])

#this goes through and plots each event using the towers
for tower_event in tEvents:
  num += 1
  grid_before = Grid(pixel_resolution=0.01, domain=domain)
  grid_before.add_tower_event(tower_event)
  grid_before.save(title='Event %d, cell resolution=0.01' % num, filename='output_tower_%d.png' % num)

  grid_after = Grid(pixel_resolution=0.01, recon_algo = 'gaussian', domain=domain)
  grid_after.add_event(tower_event.get_event())
  grid_after.save(title='Event %d, cell resolution=0.01' % num, filename='output_event_%d.png' % num)
