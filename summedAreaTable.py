import ROOT
import root_numpy as rnp
import numpy as np
import pylab as pl

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
resolution = np.array([0.1, 0.1])
#define jet radius in each direction [phi, eta]
jetRadius = np.array([1.0,1.0])

#calculate the centroidArea
# we add radius to both sides and include '1' for the center
centroidArea = np.multiply.reduce(2*jetRadius/resolution + np.array([1,1]))

#define thresholds in GeV
#  triggerable goes into our algorithm
triggerableThresh = 25.
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
  return
  print "\t"*level, 'triggerable momentum: %0.4f GeV' % (jetPt/1000.)
  print "\t"*level, 'phi: %0.4f' % (jetPhi)
  print "\t"*level, 'eta: %0.4f' % (jetEta)
  print "\t"*level, 'centroid:', tuple(phieta2cell(np.array([jetPhi,jetEta])))


def gridAddJet(grid, jetPt, jetPhi, jetEta):
  global jetRadius, resolution
  # jetCoords in cell x cell format (not phi-eta)
  jetCoords = phieta2cell(np.array([jetPhi,jetEta]))
  # amount of energy per cell inside the centroid
  jetdPt = jetPt/centroidArea
  # get tuples of the centroid's bounding box / circle
  centroidCoords = centroidMesh(jetCoords)
  #print "\t"*1, 'adding jet:'
  printCentroidDetails(jetPt, jetPhi, jetEta)
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

trigThresh = np.array([])
efficiency = np.array([])

for triggerableThresh in range(0,510,10):
  numTriggerableJets = 0
  numTriggeredJets = 0
  for event in events[[jetPt, jetPhi, jetEta]]:
    grid = gridInitialize()
    e_jetPt, e_jetPhi, e_jetEta = event
    numJets = min(e_jetPt.size, e_jetPhi.size, e_jetEta.size)
    numLargeJets = 0
    #print 'adding event with %d jets' % numJets
    for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt, e_jetPhi, e_jetEta):
      #only want jets with jetPt > 200 GeV (recorded in MeV)
      if j_jetPt/1000. < triggerableThresh:
        #print "\t"*1, 'jet is too small to be triggerable'
        printCentroidDetails(j_jetPt, j_jetPhi, j_jetEta)
        continue
      numTriggerableJets += 1
      numLargeJets += 1
      gFEXpT = gridAddJet(grid, j_jetPt, j_jetPhi, j_jetEta)
      #print "\t"*2, 'triggered momentum: %0.4f GeV' % (gFEXpT/1000.)
      if gFEXpT > triggeredThresh: numTriggeredJets += 1
    #print '-- added event with %d jets, %d triggerable' % (numJets, numLargeJets)
  print triggerableThresh, numTriggeredJets*1.0/numTriggerableJets
  trigThresh = np.append(trigThresh, triggerableThresh)
  efficiency = np.append(efficiency, numTriggeredJets*1.0/numTriggerableJets)

pl.figure()
pl.plot(trigThresh, efficiency)
pl.xlabel('$\mathrm{p}_{\mathrm{T}}^{\mathrm{jet}}$ [GeV]')
pl.ylabel('Efficiency')
pl.savefig('out.png')
