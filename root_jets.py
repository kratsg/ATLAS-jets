''' Class definitions for dealing with ATLAS jets'''

#ROOT is needed to deal with rootfiles
import ROOT

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
               pixel_resolution      = 0.20):
    '''Generates a grid of zeros based on size of domain and resolution'''
    """
        - domain: [ (phi_min, phi_max) , (eta_min, eta_max) ]
        - resolution
          - pixel: how big a single pixel is in eta-phi
    """
    self.domain                = np.array(domain).astype(float)
    self.pixel_resolution      = np.float(pixel_resolution)
    self.grid = np.zeros(self.phieta2pixel(self.__domain[:,1])).astype(float)

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

  def add_event(self, event):
    for jet in event:
      self.add_jet(jet)

  def add_jet(self, jet):
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
    return gaussian2D(pixel_jetcoord, pixel_radius, jet_energy, pixel_coord)

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
    pixel_radius = jet.radius/self.resolution
    # what we define as the jet energy for `self.__jetdPt`
    jet_energy = jet.pT
    # always start with a square mesh
    i,j = self.__square_mesh(pixel_jetcoord, pixel_radius)
    mesh = ([tuple(pixel_coord), self.__jetdPt(pixel_jetcoord, pixel_radius, jet_energy, pixel_coord)] for pixel_coord in np.transpose([i.reshape(-1), j.reshape(-1)]).astype(int) if self.boundary_conditions(coord)&circularRegion(pixel_jetcoord, pixel_radius, pixel_coord) )
    return mesh

  def __square_mesh(self, center, radius):
    '''generates a square meshgrid of points for center and sidelength 2*radius'''
    return np.meshgrid( np.arange(center[0] - radius, center[0]+radius+1), np.arange(center[1] - radius, center[1]+radius+1) )

  def __make_plot(self, title='Grid Plot'):
    '''Creates a figure of the current grid'''
    fig = pl.figure()
    # plot the grid
    pl.imshow(self.__grid.T, cmap = pl.cm.spectral)
    # x-axis is phi, y-axis is eta
    xticks_loc = pl.axes().xaxis().get_majorticklocs()
    yticks_loc = pl.axes().yaxis().get_majorticklocs()
    # make labels
    pl.xlabel('$\phi$')
    pl.ylabel('$\eta$')
    pl.title(title)
    # transform labels from pixel coords to phi-eta coords
    xticks_label = xticks_loc * self.__pixel_resolution + self.__domain[0,0]
    yticks_label = yticks_loc * self.__pixel_resolution + self.__domain[1,0]
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
    

class Jet:
  def __init__(self,\
               inputThresh    = 200.,\
               triggerThresh  = 150.,\
               E              = 0.0,\
               pT             = 0.0,\
               m              = 0.0,\
               eta            = 0.0,\
               phi            = 0.0,\
               radius         = 0.0,\
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
      radius             : radius of jet
      trigger_energy  : amount of energy actually recorded on the grid
    """
    self.inputThresh   = np.float(inputThresh)
    self.triggerThresh = np.float(triggerThresh)
    self.E             = np.float(E)
    self.pT            = np.float(pT)
    self.m             = np.float(m)
    self.coord         = (np.float(phi), np.float(eta))
    self.radius        = np.float(radius)

  def __str__(self):
    return "Jet located at (%0.4f,%0.4f) with E = %0.2f (GeV), pT = %0.2f (GeV), and m = %0.2f (GeV)"

  def triggerable(self):
    pass

  def triggered(self):
    pass

