from root_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/gFEXSlim_ttbar_zprime3000.root'

events = Events(filename=filename)
events.load()

grid = Grid(pixel_resolution=0.2)
for event in events:
  grid.add_event(event)

# at this point, all of the events have been processed (including jets)
analysis = Analysis(events)
analysis.Efficiency()
