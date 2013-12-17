from root_jets import *
filename = '/Users/kratsg/Dropbox/UChicagoSchool/DavidMiller/Data/gFEXSlim_ttbar_zprime3000.root'

events = Events(filename=filename)
events.load()

print events # Events object
print events[0] # Event object
print events[0][0] # Jet object

for event in events:
  print event
  break

for jet in events[0]:
  print jet
  break

grid = Grid()
print grid

grid.add_event(events[0])
grid.show()

grid = Grid(pixel_resolution=0.02)
for jet in events[0][0:2]:
  grid.add_jet(jet)
grid.show()
