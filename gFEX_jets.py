''' Class definitions for dealing with ATLAS jets using towers'''

#ROOT is needed to deal with rootfiles
import ROOT

#root_numpy is needed to read the rootfile
import root_numpy as rnp

# numpy and matplotlib (mpl) are used for computing and plotting 
import numpy as np
import matplotlib.pyplot as pl

def erf(x):
  '''returns the error function evaluated at x'''
  # does not necessarily need SciPy working, but less accurate without
  try:
    import scipy.special
    erf_loaded = True
  except ImportError:
    erf_loaded = False
  # use SciPy if installed
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

def circularRegion(pixel_jetcoord, pixel_radius, pixel_coord):
  '''determines if `pixel_coord` is within `pixel_radius` of `pixel_jetcoord`'''
  diff = pixel_jetcoord - pixel_coord
  distance = np.sqrt(diff[0]**2. + diff[1]**2.)
  if distance <= pixel_radius:
    return True
  else:
    return False

# all coordinates are in (phi, eta) pairs always
#     and are converted to (pixel_x, pixel_y) at runtime
class Grid:
  def __init__(self,\
               domain                = np.array([[-3.2,3.2],[0.0,3.2]]),\
               pixel_resolution      = 0.20,\
               recon_algo            = 'uniform'):
    '''Generates a grid of zeros based on size of domain and resolution'''
    """
        - domain: [ (phi_min, phi_max) , (eta_min, eta_max) ]
        - resolution
          - pixel: how big a single pixel is in eta-phi
        - recon_algo: which jet reconstruction algorithm to use
          - uniform  (implemented)
          - gaussian (implemented)
          - dipole   (not implemented)
          - j3       (not implemented)
    """
    self.domain                = np.array(domain).astype(float)
    self.pixel_resolution      = np.float(pixel_resolution)
    self.grid = np.zeros(self.phieta2pixel(self.domain[:,1])).astype(float)
    self.towers = []
    valid_algorithms = ['uniform', 'gaussian', 'dipole', 'j3']
    if recon_algo not in valid_algorithms:
      raise ValueError('%s is not a valid algorithm. Choose from %s' % (recon_algo, valid_algorithms))
    self.recon_algo = recon_algo

  def phieta2pixel(self, phieta_coord):
    '''Converts (phi,eta) -> (pixel_x, pixel_y)'''
    pixel_coord = (phieta_coord - self.domain[:,0])/self.pixel_resolution
    return np.round(pixel_coord).astype(int)

  def pixel2phieta(self, pixel_coord):
    '''Converts (pixel_x, pixel_y) -> rounded(phi, eta)'''
    return self.pixel_resolution * pixel_coord + self.domain[:,0]

  def boundary_conditions(self, pixel_coord):
    '''Checks if pixel_coord is outside of Grid'''
    # at the moment, this boundary condition defined for eta
    #   because phi is periodic
    if 0 <= pixel_coord[1] < self.grid.shape[1]:
      return True
    else:
      return False

  def add_tower_event(self, tower_event):
    if not isinstance(tower_event, TowerEvent):
      raise TypeError("You must use a TowerEvent object! You gave a %s object" % tower_event.__class__.__name__)
    for tower in tower_event:
      self.add_tower(tower)

  def add_event(self, event):
    if not isinstance(event, Event):
      raise TypeError("You must use an Event object! You gave a %s object" % event.__class__.__name__)
    for jet in event:
      self.add_jet(jet)

  def add_tower(self, tower):
    if not isinstance(tower, Tower):
      raise TypeError("You must use a Tower object! You gave us a %s object" % tower.__class__.__name__)
    '''add a single `tower` to the current grid'''
    minX, minY = self.phieta2pixel((tower.phiMin, tower.etaMin))
    maxX, maxY = self.phieta2pixel((tower.phiMax, tower.etaMax))
    print tower
    print minX, maxX, minY, maxY
    tower_mesh_coords = self.__rectangular_mesh( minX, maxX, minY, maxY )
    uniform_towerdE = tower.E/(tower_mesh_coords.size/2.)
    tower_mesh = ([tuple(pixel_coord), uniform_towerdE] for pixel_coord in tower_mesh_coords if self.boundary_conditions(pixel_coord) )
    for pixel_coord, fractional_energy in tower_mesh:
      try:
        if pixel_coord[0] >= self.grid.shape[0]:
          pixel_coord = (pixel_coord[0] - self.grid.shape[0], pixel_coord[1])
        self.grid[pixel_coord] += fractional_energy
      except IndexError:
        print "\t"*2, tower
        print "\t"*2, '-- tower pixel_coord could not be added:', pixel_coord

  def add_jet(self, jet):
    if not isinstance(jet, Jet):
      raise TypeError("You must use a Jet object! You gave us a %s object" % jet.__class__.__name__)
    '''add a single `jet` to the current grid'''
    for pixel_coord, fractional_energy in self.__generate_mesh(jet):
      try:
        if pixel_coord[0] >= self.grid.shape[0]:
          pixel_coord = (pixel_coord[0] - self.grid.shape[0], pixel_coord[1])
        self.grid[pixel_coord] += fractional_energy
        jet.trigger_energy += fractional_energy
      except IndexError:
        # we should NEVER see this gorram error
        # -- the reason is that we filter out all inappropriate eta coordinates
        #        and then wrap around in the phi coordinates
        print "\t"*2, jet
        print "\t"*2, '-- jet pixel_coord could not be added:', pixel_coord

  def __jetdPt(self, pixel_jetcoord, pixel_radius, jet_energy, pixel_coord):
    '''return the fractional energy at `pixel_coord` of a `jet_energy` GeV jet centered at `pixel_jetcoord` of radius `pixel_radius`'''
    if self.recon_algo == 'uniform':
      # uniform energy is calculated in self.__generate_mesh due to efficiency concerns
      raise Exception('should not be calling this function when self.recon_algo == \'uniform\'')
      return false
    elif self.recon_algo == 'gaussian':
      return self.__gaussian2D(pixel_jetcoord, pixel_radius, jet_energy, pixel_coord)

  def __gaussian2D(self, pixel_mu, pixel_radius, amplitude, pixel_coord):
    '''return the 2D gaussian(mu, sigma) evaluated at coord'''
    # normalization factor, 500 GeV outputs 476.275
    # default is 2 * pi * sx * sy * amplitude
    normalization = 2. * np.pi * pixel_radius**2. * erf( 0.92 * (2.**-0.5) )**2.
    exponential = np.exp(-( (pixel_mu[0] - pixel_coord[0])**2./(2. * (pixel_radius**2.)) + (pixel_mu[1] - pixel_coord[1])**2./(2.*(pixel_radius**2.)) ))
    return amplitude*exponential/normalization

  def __generate_mesh(self, jet):
    '''return the 2D mesh generator of `pixel_coords` for a `jet` to add to grid'''
    # convert to pixel coordinates and deal with grid
    pixel_jetcoord = self.phieta2pixel(jet.coord)
    pixel_radius = jet.radius/self.pixel_resolution
    # what we define as the jet energy for `self.__jetdPt`
    jet_energy = jet.pT
    # always start with a square mesh
    square_mesh_coords = self.__square_mesh(pixel_jetcoord, pixel_radius)
    if self.recon_algo == 'uniform':
      uniform_jetdPt = jet_energy/(square_mesh_coords.size/2.)
      mesh = ([tuple(pixel_coord), uniform_jetdPt] for pixel_coord in square_mesh_coords if self.boundary_conditions(pixel_coord) )
    elif self.recon_algo == 'gaussian':
      mesh = ([tuple(pixel_coord), self.__jetdPt(pixel_jetcoord, pixel_radius, jet_energy, pixel_coord)] for pixel_coord in square_mesh_coords if self.boundary_conditions(pixel_coord)&circularRegion(pixel_jetcoord, pixel_radius, pixel_coord) )
    return mesh

  def __rectangular_mesh(self, minX, maxX, minY, maxY):
    i,j = np.meshgrid( np.arange( minX, maxX), np.arange( minY, maxY) )
    return np.transpose([i.reshape(-1), j.reshape(-1)]).astype(int)

  def __square_mesh(self, center, radius):
    '''generates a square meshgrid of points for center and sidelength 2*radius'''
    return self.__rectangular_mesh(center[0] - radius, center[0] + radius + 1, center[1] - radius, center[1] + radius + 1)

  def __make_plot(self, title='Grid Plot'):
    '''Creates a figure of the current grid'''
    fig = pl.figure()
    # plot the grid
    pl.imshow(self.grid.T, cmap = pl.cm.spectral)
    # x-axis is phi, y-axis is eta
    xticks_loc = pl.axes().xaxis.get_majorticklocs()
    yticks_loc = pl.axes().yaxis.get_majorticklocs()
    # make labels
    pl.xlabel('$\phi$')
    pl.ylabel('$\eta$')
    pl.title(title)
    # transform labels from pixel coords to phi-eta coords
    xticks_label = xticks_loc * self.pixel_resolution + self.domain[0,0]
    yticks_label = yticks_loc * self.pixel_resolution + self.domain[1,0]
    pl.xticks(xticks_loc, xticks_label)
    pl.yticks(yticks_loc, yticks_label)
    # set the colorbar
    cbar = pl.colorbar(pad=0.2)
    cbar.set_label('pT (GeV)')
    return fig

  def show(self, title='Grid Plot'):
    '''Show an image of the current grid'''
    fig = self.__make_plot(title)
    fig.show()

  def save(self, title='Grid Plot', filename='output.png'):
    '''Save an image of the current grid to file'''
    fig = self.__make_plot(title)
    fig.savefig(filename)

  def __str__(self):
    return "Grid object:\n\tPhi: %s\n\tEta: %s\n\tResolution: %0.2f" % (self.domain[0], self.domain[1], self.pixel_resolution)


class Jet:
  def __init__(self,\
               inputThresh    = 200.,\
               triggerThresh  = 150.,\
               E              = 0.0,\
               pT             = 0.0,\
               m              = 0.0,\
               eta            = 0.0,\
               phi            = 0.0,\
               radius         = 1.0,\
               input_energy   = 0.0,\
               trigger_energy = 0.0):
    '''Defines a jet'''
    """
      thresholds
        - input          : in GeV, jet energy for input
        - trigger        : in GeV, jet energy for trigger
      energy             : jet energy, E
      momentum_transverse: magnitude of momentum transverse
                                 to beam, mag(p)*sin(theta)
      mass               : invariant mass of jet
      pseudo-rapidity coordinates
        - eta            : -ln( tan[theta/2] )
        - phi            : polar angle in the transverse plane
        -- theta is angle between particle momentum and the beam axis
        -- see more: http://en.wikipedia.org/wiki/Pseudorapidity
      radius             : radius of jet (eta-phi coordinates)
      input_energy       : the input energy of the jet
      trigger_energy     : amount of energy actually recorded on the grid
    """
    self.inputThresh    = np.float(inputThresh)
    self.triggerThresh  = np.float(triggerThresh)
    self.E              = np.float(E)
    self.pT             = np.float(pT)
    self.m              = np.float(m)
    self.phi            = np.float(phi)
    self.eta            = np.float(eta)
    self.coord          = (self.phi, self.eta)
    self.radius         = np.float(radius)
    self.input_energy   = np.float(pT)
    self.trigger_energy = np.float(trigger_energy)

  def __str__(self):
    return "Jet object:\n\tPhi: %0.4f\n\tEta: %0.4f\n\tE: %0.2f (GeV)\n\tpT: %0.2f (GeV)\n\tm: %0.2f (GeV)\n\tInputted: %s\n\tTriggered: %s" % (self.phi, self.eta, self.E, self.pT, self.m, self.inputted(), self.triggered())

  def inputted(self):
    return self.input_energy > self.inputThresh

  def triggered(self):
    return self.trigger_energy > self.triggerThresh

# to be grammatically correct, it should be Events' Towers
class TowerEvents:
  def __init__(self, filename = ''):
    self.filename = filename
    self.events   = []

  def __read_root_file(self):
    # read in file into a numpy record array
    self.events = rnp.root2rec(self.filename, 'minBias_14TeV_MU80_0/mytree')

  def load(self):
    self.__read_root_file()
    self.events = [TowerEvent(event=[event[i] for i in range(204,210)]) for event in self.events]

  def __iter__(self):
    # initialize to start of list
    self.iter_index = -1
    # `return self` to use `next()`
    return self

  def next(self):
    self.iter_index += 1
    if self.iter_index == len(self.events):
      raise StopIteration
    return self.events[self.iter_index]

  def __getitem__(self, index):
    return self.events[index]

  def __str__(self):
    return "TowerEvents object with %d TowerEvent objects" % len(self.events)

class TowerEvent:
  def __init__(self, event = []):
    self.towers = []
    # note that unlike David's data, it isn't a "tuple" of 215 items
    # holy mother of god, please do not blame me for the fact that
    #    I'm ignoring like 210 items in this list, we only want gTower info
    for gTowerE, gTowerNCells, gTowerEtaMin, gTowerEtaMax, gTowerPhiMin, gTowerPhiMax in zip(*event):
      self.towers.append(Tower(E=gTowerE,\
                               num_cells=gTowerNCells,\
                               etaMin=gTowerEtaMin,\
                               etaMax=gTowerEtaMax,\
                               phiMin=gTowerPhiMin,\
                               phiMax=gTowerPhiMax))

    self.E = [tower.E for tower in self.towers]
    self.phiMin = np.min([tower.phiMin for tower in self.towers])
    self.etaMin = np.min([tower.etaMin for tower in self.towers])
    self.phiMax = np.max([tower.phiMax for tower in self.towers])
    self.etaMax = np.max([tower.etaMax for tower in self.towers])

  def __iter__(self):
    # initialize to start of list
    self.iter_index = -1
    # `return self` to use `next()`
    return self

  def next(self):
    self.iter_index += 1
    if self.iter_index == len(self.towers):
      raise StopIteration
    return self.towers[self.iter_index]

  def __getitem__(self, index):
    return self.towers[index]

  def __str__(self):
    return "TowerEvents object with %d TowerEvent objects\n\tphi: (%0.4f, %0.4f)\n\teta: (%0.4f, %0.4f)" % (len(self.towers), self.phiMin, self.phiMax, self.etaMin, self.etaMax)

# to do -- fill this shit in
class Tower:
  def __init__(self,\
               E,\
               num_cells,\
               etaMin,\
               etaMax,\
               phiMin,\
               phiMax):
    self.E = E
    self.num_cells = num_cells
    self.etaMin = etaMin
    self.etaMax = etaMax
    self.phiMin = phiMin
    self.phiMax = phiMax

  def __str__(self):
    return "Tower object:\n\tE: %0.4f (GeV)\n\tnum_cells: %d\n\tphi: (%0.4f,%0.4f) \td = %0.4f\n\teta: (%0.4f, %0.4f) \td = %0.4f" % (self.E, self.num_cells, self.phiMin, self.phiMax, self.phiMax - self.phiMin, self.etaMin, self.etaMax, self.etaMax - self.etaMin)

class Events:
  def __init__(self, events = []):
    self.events = events

  def __iter__(self):
    # initialize to start of list
    self.iter_index = -1
    # `return self` to use `next()`
    return self

  def next(self):
    self.iter_index += 1
    if self.iter_index == len(self.events):
      raise StopIteration
    return self.events[self.iter_index]

  def __getitem__(self, index):
    return self.events[index]

  def __str__(self):
    return "Events object with %d Event objects" % len(self.events)
 
class Event:
  def __init__(self, event = []):
    self.jets = []
    # format generally comes as a tuple of 10 lists, each list
    #    is filled by that property for all jets like so
    #  ( [ jetE_0, jetE_1, jetE_2], [ jetPt_0, jetPt_1, jetPt_2 ], ...)
    for jetE, jetPt, jetM, jetEta, jetPhi, _, _, _, _, _ in zip(*event):
      self.jets.append(Jet(E=jetE,\
                           pT=jetPt,\
                           m=jetM,\
                           eta=jetEta,\
                           phi=jetPhi))

  def __iter__(self):
    # initialize to start of list
    self.iter_index = -1
    # `return self` to use `next()`
    return self

  def next(self):
    self.iter_index += 1
    if self.iter_index == len(self.jets):
      raise StopIteration
    return self.jets[self.iter_index]

  def __getitem__(self, index):
    return self.jets[index]

  def __str__(self):
    return "Event object with %d Jet objects" % len(self.jets)

class Analysis:
  def __init__(self, events = [], num_bins = 50):
    self.events = events
    self.num_bins = 50

  def Efficiency(self):
    input_jets = np.array([jet.pT for event in self.events for jet in [event[0], event[1]] if jet.inputted()])
    trigger_jets = np.array([jet.pT for event in self.events for jet in [event[0], event[1]] if jet.triggered()])

    bin_range = (input_jets.min(), input_jets.max())
    histogram_input = np.histogram(input_jets, range=bin_range, bins=self.num_bins)
    histogram_trigger = np.histogram(trigger_jets, range=bin_range, bins=self.num_bins)
    nonzero_bins = np.where(histogram_input[0] != 0)
    efficiency = np.true_divide(histogram_trigger[0][nonzero_bins], histogram_input[0][nonzero_bins])

    pl.figure()
    pl.scatter(histogram_input[1][nonzero_bins], efficiency)
    pl.xlabel('$\mathrm{p}_{\mathrm{T}}^{\mathrm{jet}}$ [GeV]')
    pl.ylabel('Efficiency')
    pl.title('Turn-on Plot')
    pl.show()
