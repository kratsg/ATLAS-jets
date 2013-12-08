import ROOT
import root_numpy as rnp
import numpy as np
import pylab as pl
import matplotlib

#absolute filenames make sense, easier to find what file we talking about
filename = '/Users/giordon/Dropbox/UChicagoSchool/DavidMiller/Data/gFEXSlim_ttbar_zprime3000.root'
#load up the root file into a numpy.darray
events = rnp.root2array(filename, 'jets')
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
resolution = np.array([0.2, 0.2])
#define jet radius in each direction [phi, eta]
jetRadius = np.array([1.0,1.0])
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

#generate a meshgrid for a given jet centroid
def centroidMesh(jetCoord):
  global jetRadius, resolution
  i,j = np.meshgrid(np.arange(jetCoord[0]-jetRadius[0]/resolution[0], jetCoord[0]+jetRadius[0]/resolution[0]+1), np.arange(jetCoord[1]-jetRadius[1]/resolution[1], jetCoord[1]+jetRadius[1]/resolution[1]+1))
  return map(lambda x: tuple(x), np.transpose([i.reshape(-1), j.reshape(-1)]).astype(int))

#define boolean function for allowed coordinates (wrap-around, etc)
def boundaryConditions(grid,coord):
  #coord = (phi, eta) in cell coordinates
  if 0 <= coord[1] < grid.shape[1]:
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
  # amount of energy per cell inside the centroid
  jetdPt = jetPt/(1000.*centroidArea)
  # get tuples of the centroid's bounding box / circle
  centroidCoords = centroidMesh(jetCoords)
  #print "\t"*1, 'adding jet:'
  #printCentroidDetails(jetPt, jetPhi, jetEta)
  # number of jet points
  numJetCoords = len(centroidCoords)
  # number of recorded points
  numJetCoordsActual = 0
  # now to loop over all coordinates in the centroid that fit
  for coord in filter(lambda x: boundaryConditions(grid,x), centroidCoords):
    # handle periodicity in here! more specifically, if phi is too large, it doesn't wrap around
    try:
      numJetCoordsActual += 1
      if coord[0] >= grid.shape[0]:
        coord = (coord[0] - grid.shape[0], coord[1])
      grid[coord] += jetdPt
    except IndexError:
      # we should NEVER see this gorram error
      # -- the reason is that we filter out all inappropriate eta coordinates
      #        and then wrap around in the phi coordinates
      print "\t"*2, '-- jet coord could not be added:', coord
  #print "\t"*3, "%d/%d points added ~ %0.2f%%" % (numJetCoordsActual, numJetCoords, numJetCoordsActual*100.0/numJetCoords)
  return numJetCoordsActual*jetPt/numJetCoords

def gridPlotEvent(grid, triggerableJets):
  pl.figure()
  print "\t"*1,"making plot for event"
  #adding fake elements to fix legend
  fake = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none', visible=False)
  for jetPt, jetPhi, jetEta in triggerableJets:
    print "\t"*2, "(%0.2f, %0.2f) -- %0.3f" % (jetPhi,jetEta,jetPt/1000.)
  #transpose is necessary! want x = eta, y = phi
  pl.imshow(grid.T, cmap = pl.cm.spectral)
  #x is phi
  xticks_loc = pl.axes().xaxis.get_majorticklocs()
  #y is eta
  yticks_loc = pl.axes().yaxis.get_majorticklocs()
  #top two jet energies
  topTwoEnergies = triggerableJets[:,0].argsort()[::-1][:2]
  topEnergy, nextTopEnergy = np.round(triggerableJets[:,0][topTwoEnergies]/1000.)
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
  pl.show()


def computeEfficiency():
  global events, triggerableThresh
  triggerableJets = []
  triggeredJets = []
  for event in events[[jetPt, jetPhi, jetEta]]:
    #triggerableJets = []
    grid = gridInitialize()
    e_jetPt, e_jetPhi, e_jetEta = event
    numJets = e_jetPt.size
    #we only want to compute efficiency for top two jets
    #topTwoJets = e_jetPt.argsort()[::-1][:2]
    #print 'adding event with %d jets' % numJets
    for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt, e_jetPhi, e_jetEta):
    #for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt[topTwoJets], e_jetPhi, e_jetEta):
      #only want jets with jetPt > 200 GeV (recorded in MeV)
      if j_jetPt/1000. < triggerableThresh:
        continue
      triggerableJets.append([j_jetPt, j_jetPhi, j_jetEta])
      gFEXpT = gridAddJet(grid, j_jetPt, j_jetPhi, j_jetEta)
      if gFEXpT/1000. > triggeredThresh:
        triggeredJets.append([j_jetPt, j_jetPhi, j_jetEta])
    #if numJets > 3 or e_jetPt.max()/1000. >= 600.:
    #  gridPlotEvent(grid, np.array(triggerableJets))
  triggerableJets = np.array(triggerableJets)
  triggeredJets = np.array(triggeredJets)
  binRange = (triggerableJets[:,0].min(), triggerableJets[:,0].max())
  numBins = 50
  histTriggerable = np.histogram(triggerableJets[:,0], range=binRange, bins=numBins)
  histTriggered = np.histogram(triggeredJets[:,0], range=binRange, bins=numBins)
  nonZerobins = np.where(histTriggerable[0] != 0)
  efficiency = np.true_divide(histTriggered[0][nonZerobins], histTriggerable[0][nonZerobins])
  pl.figure()
  pl.scatter(histTriggerable[1][nonZerobins]/1000., efficiency)
  pl.xlabel('$\mathrm{p}_{\mathrm{T}}^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Efficiency')
  pl.title('%d bins, Triggerable = %d GeV, Triggered = %d GeV' % (numBins, triggerableThresh, triggeredThresh))
  pl.savefig('efficiency_jets_%d_%d.png' % (triggerableThresh, triggeredThresh))
  #now compute the efficiency inside the gFEX by applying boundary conditions to filter out centroid locations
  gFEX_triggerableJets = np.array(filter(lambda x: boundaryConditions(grid,phieta2cell((x[1],x[2]))), triggerableJets))
  gFEX_removedJets = np.array(filter(lambda x: not boundaryConditions(grid,phieta2cell((x[1],x[2]))), triggerableJets))
  gFEX_triggeredJets = np.array([jet for jet in triggeredJets if jet not in gFEX_removedJets])
  gFEX_histTriggerable = np.histogram(gFEX_triggerableJets[:,0], range=binRange, bins=numBins)
  gFEX_histTriggered = np.histogram(gFEX_triggeredJets[:,0], range=binRange, bins=numBins)
  gFEX_nonZerobins = np.where(gFEX_histTriggerable[0] != 0)
  gFEX_efficiency = np.true_divide(gFEX_histTriggered[0][gFEX_nonZerobins], gFEX_histTriggerable[0][gFEX_nonZerobins])
  pl.figure()
  pl.scatter(gFEX_histTriggerable[1][gFEX_nonZerobins]/1000., gFEX_efficiency)
  pl.xlabel('$\mathrm{p}_{\mathrm{T}}^{\mathrm{jet}}$ [GeV]')
  pl.ylabel('Efficiency')
  pl.title('%d bins, Triggerable = %d GeV, Triggered = %d GeV' % (numBins, triggerableThresh, triggeredThresh))
  pl.savefig('efficiency_gFEX_%d_%d.png' % (triggerableThresh, triggeredThresh))

computeEfficiency()
