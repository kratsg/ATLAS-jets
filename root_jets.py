''' Class definitions for dealing with ATLAS jets'''

#ROOT is needed to deal with rootfiles
import ROOT

# numpy and matplotlib (mpl) are used for computing and plotting 
import numpy as np
import matplotlib.pyplot as pl

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
    self.__domain                = np.array(domain).astype(float)
    self.__pixel_resolution      = np.float(pixel_resolution)
    self.__grid = np.zeros(self.phieta2pixel(self.__domain[:,1])).astype(float)

  def phieta2pixel(self, phieta_coord):
    '''Converts (phi,eta) -> (pixel_x, pixel_y)'''
    pixel_coord = (phieta_coord - self.__domain[:,0])/self.__pixel_resolution
    return np.round(pixel_coord).astype(int)

  def pixel2phieta(self, pixel_coord):
    '''Converts (pixel_x, pixel_y) -> rounded(phi, eta)'''
    return self.__pixel_resolution * pixel_coord + self.__domain[:,0]

  def boundary_conditions(self, pixel_coord):
    '''Checks if pixel_coord is outside of Grid'''
    # at the moment, this boundary condition defined for eta
    #   because phi is periodic
    if 0 <= pixel_coord[1] < grid.shape[1]:
      return True
    else:
      return False

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
    
