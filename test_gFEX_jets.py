from gFEX_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim_TTbar_14TeV_MU80.root'
directory = 'TTbar_14TeV_MU80/'
tEvents = TowerEvents(filename=filename, directory=directory, seed_filter=SeedFilter(ETthresh=20.0))
tEvents.load()

print 'Loaded and built all objects for gTowers data'
num = 0

domain = np.array([[-3.2,3.2],[-4.9,4.9]])

#this goes through and plots each event using the towers
for tower_event in tEvents:
  num += 1
  grid_before = Grid(pixel_resolution=0.01, domain=domain)
  grid_before.add_tower_event(tower_event)
  grid_before.save(title='Event %d, cell resolution=0.01' % num, filename='output_tower_%d.png' % num, colzLabel = '$E_T^{\mathrm{tower}}$')

  grid_after = Grid(pixel_resolution=0.01, recon_algo = 'gaussian', domain=domain)
  grid_after.add_event(tower_event.get_event())
  grid_after.save(title='Event %d, cell resolution=0.01' % num, filename='output_event_%d.png' % num, colzLabel = '$p_T^{\mathrm{jet}}$')

i = 0
bins = np.arange(0,15,0.25)
cumul_sum = np.zeros(len(bins)-1)
for tower_event in tEvents:
  #this goes ahead and makes histograms of the ET distro of the towers for each event
  pl.figure()
  ETs = [tower.E/np.cosh(tower.eta) for tower in tower_event.towers]
  pl.hist(ETs, bins=bins)
  pl.xlabel('$E_T^{\mathrm{tower}}$ [GeV]')
  pl.ylabel('Number of gTowers')
  i += 1
  pl.title('Histogram of $E_T^{\mathrm{tower}}$ for Event #%d, %d towers' % (i, len(tower_event.towers)))
  pl.xlim(0,15)
  pl.ylim(0,400)
  pl.savefig('event_%d_histogram.png' % i)
  pl.close()
  #now we want the histogram of towerThresholds, so we need to count # towers above a certain threshold of ET
  hist, bin_edges = np.histogram(ETs, bins=bins)
  cumul_sum += np.sum(hist) - np.cumsum(hist)

pl.figure()
pl.xlabel('$E_T^{\mathrm{threshold}}$ [GeV]')
pl.ylabel('Number of gTowers')
pl.title('Number of gTowers above $E_T^{\mathrm{threshold}}$')
pl.bar(bin_edges[:-1], cumul_sum, width=0.25)
pl.xlim(0,15)
pl.savefig('events_threshold_histogram.png')
pl.close()
