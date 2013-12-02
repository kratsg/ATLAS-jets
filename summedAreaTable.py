import ROOT
import root_numpy as rnp
import numpy as np

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
  print "\t"*level, 'transverse momentum: %0.4f GeV' % (jetPt/1000.)
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
  print "\t"*1, 'adding jet:'
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
  print "\t"*2, '-- jet centroid added'
  print "\t"*3, "%d/%d points added ~ %0.2f%%" % (numJetCoordsActual, numJetCoords, numJetCoordsActual*100.0/numJetCoords)

for event in events[[jetPt, jetPhi, jetEta]]:
  grid = gridInitialize()
  e_jetPt, e_jetPhi, e_jetEta = event
  numJets = min(e_jetPt.size, e_jetPhi.size, e_jetEta.size)
  numSmallJets = 0
  print 'adding event with %d jets' % numJets
  for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt, e_jetPhi, e_jetEta):
    #only want jets with jetPt > 200 GeV (recorded in MeV)
    if j_jetPt/1000. < 200.:
      print "\t"*1, 'jet is too small to trigger'
      printCentroidDetails(j_jetPt, j_jetPhi, j_jetEta)
      numSmallJets += 1
      continue
    gridAddJet(grid, j_jetPt, j_jetPhi, j_jetEta)
  print '-- added event with %d jets, %d jets did not trigger, %d added to grid' % (numJets, numSmallJets, numJets - numSmallJets)
  break


