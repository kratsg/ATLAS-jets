Repo containing code for analyzing jets from ATLAS

Current Features
====================

- adds jets to a grid that satisfy jetPt > 200. GeV and are within the grid's boundaries of `[eta,phi] = [(0,+3.2), (-3.2,3.2)]`
- reports portions of the centroid that are outside said boundary (with the specific coordinates)
- reports metrics about centroid coordinates that did not get added to grid
- reports how many jets are in an event
- reports how many jets for an event did not trigger (above threshold of 200 GeV)
- reports how many jets for an event added to grid
- reports how many points for a jet are added to grid
- wraps points around in \phi but not in \eta

### To Do

