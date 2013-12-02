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


Example Output (for event)
===========================
```
adding event with 8 jets
  adding jet:
    transverse momentum: 1398.3260 GeV
    phi: 0.5748
    eta: -0.3457
    centroid: (19, -2)
    -- jet centroid added
      44/121 points added ~ 36.36%
  adding jet:
    transverse momentum: 1226.3349 GeV
    phi: -2.3851
    eta: -0.3746
    centroid: (4, -2)
    -- jet centroid added
      44/121 points added ~ 36.36%
  adding jet:
    transverse momentum: 256.8598 GeV
    phi: 2.5848
    eta: -0.0515
    centroid: (29, 0)
    -- jet centroid added
      66/121 points added ~ 54.55%
  jet is too small to trigger
    transverse momentum: 32.1432 GeV
    phi: -1.1973
    eta: 0.3768
    centroid: (10, 2)
  jet is too small to trigger
    transverse momentum: 20.3603 GeV
    phi: 3.1049
    eta: -1.9862
    centroid: (32, -10)
  jet is too small to trigger
    transverse momentum: 18.5656 GeV
    phi: 1.2787
    eta: 2.3406
    centroid: (22, 12)
  jet is too small to trigger
    transverse momentum: 13.6911 GeV
    phi: -1.2916
    eta: -1.9791
    centroid: (10, -10)
  jet is too small to trigger
    transverse momentum: 11.7000 GeV
    phi: 0.8998
    eta: -4.0269
    centroid: (20, -20)
-- added event with 8 jets, 5 jets did not trigger, 3 added to grid

```
