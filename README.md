# DirectPhotonAnalysis

## EmcalEtViewer
### DisplayLegoPlot.cc (Fun4All module)
This subsysreco module dumps the E_T values for all EMCal towers in a ttree named `cemc_et`.

The E_T value is calculated as `ET = tower->get_energy() / cosh(eta_from_vertex)`.
The eta is evaluated

### DrawCemcLego
This ROOT macro reads a ttree and saves an image file (.png) of Lego plot event by event.

Usage:
```
[]$ root 'DrawCemcLego.C("/sphenix/user/ryotaro/DisplayLegoPlot/output/0/0-00000.root")'
```
