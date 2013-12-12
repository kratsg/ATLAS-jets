import ROOT
import numpy as np
import pylab as pl
import matplotlib

events = ''
def read_file():
  global events
  try:
    import root_numpy as rnp
    # absolute filenames make sense, easier to find what file we talking about
    filename = '/Users/giordon/Dropbox/UChicagoSchool/DavidMiller/Data/gFEXSlim_ttbar_zprime3000.root'
    # load up the root file into a numpy.recarray
    events = rnp.root2rec(filename, 'jets')
    print "read file using root_numpy"
  except:
    import pickle
    filename = open('jetsData.pkl','rb')
    events = pickle.load(filename)
    print "read file using pickle"

read_file()
#map pointers to column names within the file
jetE,jetPt,jetM,jetEta,jetPhi,_,_,_,_,_ = events.dtype.names
#jetE: jet energy
#jetPt: jet transverse momentum
#jetM: jet mass
#jetEta: jet eta coordinate
#jetPhi: jet phi coordinate

#at this point, we have a structure like
'''
events: [
  (jetE, jetPt, jetM, jetEta, ...),
  (jetE, jetPt, jetM, jetEta, ...),
  (jetE, jetPt, jetM, jetEta, ...)
]
and each tuple contains a list for each jet. Not the most convenient data access

For each event, we want to generate an eti-phi grid of jets
  implementing a variable radius corresponding to the center
  of each jet location, and then apply a summed area table
  for that event's grid to determine if we would throw away
  that event or not!

'''
#phi == x, eta = y
#want domain of grid, not of the total window view
domain = np.array([[-3.2,3.2], [0.0, 3.2]])
#resolution is basically gridSize for right now
resolution = np.array([0.05, 0.05])
#define jet radius in each direction [phi, eta]
jetRadius = np.array([0.85,0.85])
#calculate the centroidArea
# we add radius to both sides and include '1' for the center
centroidArea = np.multiply.reduce(2*jetRadius/resolution + np.array([1,1]))
#define thresholds in GeV
#  triggerable goes into our algorithm
triggerableThresh = 150.
#  triggered is what would be triggered by our algorithm
triggeredThresh = 200.

#converts (phi,eta) -> (cellPhi, cellEta) coordinates
def phieta2cell(coord):
  return np.round((coord - domain[:,0])/resolution).astype(int)

#converts (cellPhi, cellEta) -> (phi,eta) coordinate
def cell2phieta(coord):
  return resolution*coord + domain[:,0]

#initialize the viewing window grid
def gridInitialize():
  global domain
  #note that phieta2cell(domain[:,0]) = (0,0) by definition
  return np.zeros(phieta2cell(domain[:,1])).astype(float)

# work around since David doesn't have SciPy working
# will try to use scipy.special.erf(x), but defaults to
#       custom function (not as accurate)
def erf(x):
  try: 
    import scipy.special
    erf_loaded = True 
  except ImportError: 
    erf_loaded = False 

  if erf_loaded:
    return scipy.special.erf(x)
  else:
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = np.fabs(x)
    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911
    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

# returns a 2D gaussian centered at mu evaluated at coord (must be in cell grid)
def gaussian2D(mu, sigma, amplitude, coord):
  # normalize gaussian so that the integral is just the amplitude
  #    Note: normalization is 2 pi sx sy * amplitude

  # normalized such that testScript outputs 476.275
  A = amplitude/(2*np.pi*sigma[0]*sigma[1]*erf(0.92*(2.**-0.5))**2.)
  #A = amplitude/(2*np.pi*sigma[0]*sigma[1]*scipy.special.erf(2.**-0.5)**2.)
  exponential = np.exp(-( (mu[0] - coord[0])**2./(2.*(sigma[0]**2.)) + (mu[1] - coord[1])**2./(2.*(sigma[1]**2.))  ))
  return A*exponential

#returns the fractional energy at the given coordinates for a jet
def jetdPt(jetCoord, jetEnergy, coord):
  global jetRadius, resolution
  return gaussian2D(jetCoord, jetRadius/resolution, jetEnergy, coord)

#generate a meshgrid for a given jet centroid
# jetCoord in cell coordinates
def centroidMesh(grid, jetCoord, jetPt):
  global jetRadius, resolution
  i,j = np.meshgrid(np.arange(jetCoord[0]-jetRadius[0]/resolution[0], jetCoord[0]+jetRadius[0]/resolution[0]+1), np.arange(jetCoord[1]-jetRadius[1]/resolution[1], jetCoord[1]+jetRadius[1]/resolution[1]+1))
  mesh = map(lambda x: [tuple(x), jetdPt(jetCoord, jetPt, x)], np.transpose([i.reshape(-1), j.reshape(-1)]).astype(int))
  mesh = filter(lambda x: boundaryConditions(grid,x[0])&circularRegion(jetRadius/resolution, jetCoord, x[0]), mesh)
  return mesh

#define boolean function for allowed coordinates (wrap-around, etc)
def boundaryConditions(grid, coord):
  #coord = (phi, eta) in cell coordinates
  if 0 <= coord[1] < grid.shape[1]:
    return True
  else:
    return False

def circularRegion(radius, jetCoord, coord):
  #coord = (phi, eta) in cell coordinates
  diff = jetCoord - coord
  distance = np.sqrt(diff[0]**2. + diff[1]**2.)
  if distance <= radius[0]:
    return True
  else:
    return False

#quickly print details about a jet
def printCentroidDetails(jetPt, jetPhi, jetEta, level=2):
  print "\t"*level, 'triggerable momentum: %0.4f GeV' % (jetPt/1000.)
  print "\t"*level, 'phi: %0.4f' % (jetPhi)
  print "\t"*level, 'eta: %0.4f' % (jetEta)
  print "\t"*level, 'centroid:', tuple(phieta2cell(np.array([jetPhi,jetEta])))


def gridAddJet(grid, jetPt, jetPhi, jetEta):
  global jetRadius, resolution
  # jetCoords in cell x cell format (not phi-eta)
  jetCoords = phieta2cell(np.array([jetPhi,jetEta]))

  # get tuples of the centroid's bounding box / circle
  # this is in the form of [(x,y), fractionalEnergy]
  centroidCoords = centroidMesh(grid, jetCoords, jetPt)
  # jet energy added to grid
  energyRecorded = 0.0
  # now to loop over all coordinates in the centroid that fit
  # coordinates are of the form [(x,y), fractional energy]
  for coord, fractionalEnergy in centroidCoords:
    # handle periodicity in here! more specifically, if phi is too large, it doesn't wrap around
    try:
      if coord[0] >= grid.shape[0]:
        coord = (coord[0] - grid.shape[0], coord[1])
      grid[coord] += fractionalEnergy
      energyRecorded += fractionalEnergy
    except IndexError:
      # we should NEVER see this gorram error
      # -- the reason is that we filter out all inappropriate eta coordinates
      #        and then wrap around in the phi coordinates
      print "\t"*2, '-- jet coord could not be added:', coord
  #print jetPt, energyRecorded, cell2phieta(jetCoords)
  return energyRecorded

def gridPlotEvent(grid, triggerableJets):
  if len(triggerableJets) < 2:
    print "no jets to see (or only one...)"
    return
  pl.figure()
  print "\t"*1,"making plot for event"
  #adding fake elements to fix legend
  fake = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none', visible=False)
  for jetPt, jetPhi, jetEta in triggerableJets:
    print "\t"*2, "(%0.2f, %0.2f) -- %0.3f" % (jetPhi,jetEta,jetPt)
  #transpose is necessary! want x = eta, y = phi
  pl.imshow(grid.T, cmap = pl.cm.spectral)
  #x is phi
  xticks_loc = pl.axes().xaxis.get_majorticklocs()
  #y is eta
  yticks_loc = pl.axes().yaxis.get_majorticklocs()
  #top two jet energies
  topTwoEnergies = triggerableJets[:,0].argsort()[::-1][:2]
  topEnergy, nextTopEnergy = np.round(triggerableJets[:,0][topTwoEnergies])
  pl.xlabel('$\phi$')
  pl.ylabel('$\eta$')
  pl.title('Event - %d triggerable jets. Leading energies (GeV): %d, %d' % (len(triggerableJets), topEnergy, nextTopEnergy))
  pl.legend([fake]*len(triggerableJets), ['(%0.2f, %0.2f)' % (jetPhi,jetEta) for jetPt, jetPhi, jetEta in triggerableJets], prop={'size':12}, loc="center", bbox_to_anchor=(0.5,-0.5), ncol=4, fancybox=True, shadow=True, title='Jets')
  #transform labels from cellGrid coords to phi-eta coords
  xticks_label = xticks_loc*resolution[0] + domain[0,0]
  yticks_label = yticks_loc*resolution[1] + domain[1,0]
  pl.xticks(xticks_loc, xticks_label)
  pl.yticks(yticks_loc, yticks_label)
  cbar = pl.colorbar(shrink=0.75, pad=.2)
  cbar.set_label('pT (GeV)')
  pl.savefig("%d-jets_top-%d-GeV_next-%d-GeV.png" % (len(triggerableJets), np.round(topEnergy), np.round(nextTopEnergy)))

def computeEfficiency():
  global events, triggerableThresh
  triggerableJets = []
  triggeredJets = []
  i = 0
  for event in events:
    i+= 1
    if i % 25 == 0:
      print "Event #%d" % i
    grid = gridInitialize()
    e_jetE,e_jetPt,e_jetM,e_jetEta,e_jetPhi,_,_,_,_,_ = event
    numJets = e_jetPt.size
    #we only want to compute efficiency for top two jets
    #topTwoJets = e_jetPt.argsort()[::-1][:2]
    #print 'adding event with %d jets' % numJets
    for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt, e_jetPhi, e_jetEta):
      #only want jets with jetPt > 200 GeV (recorded in MeV)
      if j_jetPt/1000. < triggerableThresh:
        continue
      triggerableJets.append([j_jetPt/1000., j_jetPhi, j_jetEta])
      gFEXpT = gridAddJet(grid, j_jetPt/1000., j_jetPhi, j_jetEta)
      if gFEXpT > triggeredThresh:
        triggeredJets.append([j_jetPt/1000., j_jetPhi, j_jetEta])
  triggerableJets = np.array(triggerableJets)
  triggeredJets = np.array(triggeredJets)
  binRange = (triggerableJets[:,0].min(), triggerableJets[:,0].max())
  numBins = 50
  histTriggerable = np.histogram(triggerableJets[:,0], range=binRange, bins=numBins)
  histTriggered = np.histogram(triggeredJets[:,0], range=binRange, bins=numBins)
  nonZerobins = np.where(histTriggerable[0] != 0)
  efficiency = np.true_divide(histTriggered[0][nonZerobins], histTriggerable[0][nonZerobins])
  pl.figure()
  pl.scatter(histTriggerable[1][nonZerobins], efficiency)
  pl.xlabel('$\mathrm{p}_{\mathrm{T}}^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Efficiency')
  pl.title('%d bins, Triggerable = %d GeV, Triggered = %d GeV' % (numBins, triggerableThresh, triggeredThresh))
  pl.savefig('efficiency_jets_%d_%d.png' % (triggerableThresh, triggeredThresh))
  #now compute the efficiency inside the gFEX by applying boundary conditions to filter out centroid locations

  gFEX_triggerableJets = np.array(filter(lambda x: x[2] >= 0.0, triggerableJets))
  gFEX_removedJets = np.array(filter(lambda x: x[2] < 0.0, triggerableJets))
  gFEX_triggeredJets = np.array([jet for jet in triggeredJets if jet not in gFEX_removedJets])
  gFEX_histTriggerable = np.histogram(gFEX_triggerableJets[:,0], range=binRange, bins=numBins)
  gFEX_histTriggered = np.histogram(gFEX_triggeredJets[:,0], range=binRange, bins=numBins)
  gFEX_nonZerobins = np.where(gFEX_histTriggerable[0] != 0)
  gFEX_efficiency = np.true_divide(gFEX_histTriggered[0][gFEX_nonZerobins], gFEX_histTriggerable[0][gFEX_nonZerobins])
  pl.figure()
  pl.scatter(gFEX_histTriggerable[1][gFEX_nonZerobins], gFEX_efficiency)
  pl.xlabel('$\mathrm{p}_{\mathrm{T}}^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Efficiency')
  pl.title('%d bins, Triggerable = %d GeV, Triggered = %d GeV' % (numBins, triggerableThresh, triggeredThresh))
  pl.savefig('efficiency_gFEX_%d_%d.png' % (triggerableThresh, triggeredThresh))

def showEvents(maxNum=25):
  global events, triggerableThresh
  i = 0
  for event in events:
    triggerableJets = []
    triggeredJets = []
    grid = gridInitialize()
    e_jetE,e_jetPt,e_jetM,e_jetEta,e_jetPhi,_,_,_,_,_ = event
    #we only want to compute efficiency for top two jets
    #topTwoJets = e_jetPt.argsort()[::-1][:2]
    #print 'adding event with %d jets' % numJets
    for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt, e_jetPhi, e_jetEta):
    #for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt[topTwoJets], e_jetPhi, e_jetEta):
      #only want jets with jetPt > 200 GeV (recorded in MeV)
      if j_jetPt/1000. < triggerableThresh:
        continue
      triggerableJets.append([j_jetPt/1000., j_jetPhi, j_jetEta])
      gFEXpT = gridAddJet(grid, j_jetPt/1000., j_jetPhi, j_jetEta)
      if gFEXpT > triggeredThresh:
        triggeredJets.append([j_jetPt/1000., j_jetPhi, j_jetEta])
    gridPlotEvent(grid, np.array(triggerableJets))
    print "\t Triggered Jets"
    for triggeredJet in triggeredJets:
      print "\t\t Energy:", triggeredJet[0], "\t", "Coords:", (triggeredJet[1],triggeredJet[2])
    i+=1
    if i >= maxNum:
      return True


def testScript():
  global resolution
  resolution /= 10.
  grid = gridInitialize()
  print gridAddJet(grid, 500., 0.0, 1.6)

print "Running compute efficiency at resolution of 0.2 x 0.2"
resolution *= 10.
computeEfficiency()
print "Now running showEvents(5)"
showEvents(5)
