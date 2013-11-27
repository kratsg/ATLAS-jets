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

def gridInitialize():
  global domain, resolution
  #determine size of grid
  x,y = domain/resolution
  rows, cols = np.array([ x.max() - x.min(), y.max() - y.min() ])
  return np.zeros((rows,cols))


def gridAddJet(grid, jetPt, jetPhi, jetEta):
  global jetRadius
  # jetCoords in cell x cell format (not phi-eta)
  jetCoords = np.round(np.array([jetPhi, jetEta])/resolution).astype(int)
  # amount of energy per cell inside the centroid
  jetdPt = jetPt/centroidArea
  # now to loop over all coordinates in the centroid
  #   first make a meshgrid of coordinates in the centroid
  i,j = np.meshgrid(np.arange(jetCoords[0]-jetRadius[0], jetCoords[0]+jetRadius[0]+1), np.arange(jetCoords[1]-jetRadius[1], jetCoords[1]+jetRadius[1]+1))
  # reshape to a flattened array so we get pairs like (0,0), (0,1), ...
  centroidCoords = np.transpose([i.reshape(1,-1), j.reshape(1,-1)])
  print 'adding jet'
  print "\t", jetPt, jetPhi, jetEta
  for coord in centroidCoords:
    try:
      grid[coord] = jetdPt
    except IndexError:
      print 'Jet centroid coordinates are out of range'
      print "\t\t", coord

for event in events[[jetPt, jetPhi, jetEta]]:
  grid = gridInitialize()
  e_jetPt, e_jetPhi, e_jetEta = event
  for j_jetPt, j_jetPhi, j_jetEta in zip(e_jetPt, e_jetPhi, e_jetEta):
    #only want jets with jetPt > 200 GeV (recorded in MeV)
    if j_jetPt/1000. < 200.:
      continue
    gridAddJet(grid, j_jetPt, j_jetPhi, j_jetEta)

