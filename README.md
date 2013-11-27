Repo containing code for analyzing jets from ATLAS

Current Features
====================

- adds jets to a grid that satisfy jetPt > 200. GeV and are within the grid's boundaries of `[eta,phi] = [(0,+3.2), (-3.2,3.2)]`
- reports portions of the centroid that are outside said boundary (with the specific coordinates)

### To Do

- add periodicity in \phi
- add metrics for each event such as number of jets and so forth
