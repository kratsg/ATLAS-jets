from gFEX_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/PileupSkim_TTbar_14TeV_MU80.root'
directory = 'TTbar_14TeV_MU80/'
tEvents = TowerEvents(filename=filename, directory=directory, seed_filter=SeedFilter(ETthresh=20.0))
tEvents.load()

print 'Loaded and built all objects for gTowers data'
#events = Events([tower_event.get_event() for tower_event in tEvents])

num = 0

domain = np.array([[-3.2,3.2],[-4.9,4.9]])

i = 0
num_bins = 100
cumul_sum = np.zeros(num_bins).astype(int)
# computing maxET (needed for consistent binning)
maxET = np.max([np.max([tower.E/np.cosh(tower.eta) for tower in towers]) for towers in tEvents])
maxET += 0.1 #this is so we can fix binning issues
bin_edges = np.arange(0.0, maxET + maxET/num_bins, maxET/num_bins)

for tower_event in tEvents:
  #this goes ahead and makes histograms of the ET distro of the towers for each event
  #pl.figure()
  ETs = [tower.E/np.cosh(tower.eta) for tower in tower_event.towers]
  #pl.hist(ETs, bins=bins)
  #pl.xlabel('$E_T^{\mathrm{tower}}$ [GeV]')
  #pl.ylabel('Number of gTowers')
  #i += 1
  #pl.title('Histogram of $E_T^{\mathrm{tower}}$ for Event #%d, %d towers' % (i, len(tower_event.towers)))
  #pl.savefig('event_%d_histogram.png' % i)
  #pl.close()
  #now we want the histogram of towerThresholds, so we need to count # towers above a certain threshold of ET
  hist, _ = np.histogram(ETs, bins=bin_edges)
  cumul_sum += np.cumsum(hist[::-1])[::-1]

#get width of each bar based on domain/num_bins
width=(bin_edges[-1] - bin_edges[0])/(num_bins+5.)
#normalize distribution per Michael Begel
cumul_sum = 1.0*cumul_sum/len(tEvents.events)
# plot it all
pl.figure()
pl.xlabel('$E_T^{\mathrm{threshold}}$ [GeV]')
pl.ylabel('Number of gTowers')
pl.title('Number of gTowers above $E_T^{\mathrm{threshold}}$')
pl.bar(bin_edges[:-1], cumul_sum, width=width, log=True)
#pl.xscale('log')
#pl.yscale - need to use log=True argument in pyplot.bar (see documentation)
pl.savefig('events_threshold_histogram.png')
pl.close()

print "DONE WITH HISTOGRAMS"
raise SystemExit

#this goes through and plots each event using the towers
for tower_event in tEvents:
  num += 1
  grid_before = Grid(pixel_resolution=0.01, domain=domain)
  grid_before.add_tower_event(tower_event)
  grid_before.save(title='Event %d, cell resolution=0.01' % num, filename='output_tower_%d.png' % num, colzLabel = '$E_T^{\mathrm{tower}}$')

  grid_after = Grid(pixel_resolution=0.01, recon_algo = 'gaussian', domain=domain)
  grid_after.add_event(tower_event.get_event())
  grid_after.save(title='Event %d, cell resolution=0.01' % num, filename='output_event_%d.png' % num, colzLabel = '$p_T^{\mathrm{jet}}$')


