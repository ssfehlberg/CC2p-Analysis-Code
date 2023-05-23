# CC2p-Analysis-Code
Ths repository contains a variety of code used 

## RooUnfold
This folder is a clone of the [RooUnfold repository]([url](https://gitlab.cern.ch/RooUnfold/RooUnfold)). We make use of RooUnfold's implementation of D'Agostini to unfold the data from reconstructed space to true space. Note: the files in this folder are specific to my personal laptop i.e. individuals will not be able to run the following code without a system specific copy of RooUnfold. It should be as simple as redoing the make. 

## PeLEE
This folder contains code to create a variety of different plots using the output root files created from code found in the [Event Selection Repository]([url](https://github.com/ssfehlberg/CC2p-Event-Selection)). Produced plots include 1) efficiency distributions 2) X,Y, and Z coordinates of the reconstructed vertex and 3) the selected event distributions. There are also commented out blocks of code used to generate a variety of PFP plots and plots of cuts (such as the PID and track score value). 

```
root -b analysis.C
analysis s
s.main()
```
Users will then be given two different prompts. The first will indicate if you wish to look at the distributions using the nomial binning (0 = pelee) or with the optimized cross section binning (1 = pelee with xsec binning). The second prompt will ask which run you wish to process. 

### A few notes:
- If you have not yet run the systematic code, see below first.
- You can find all the various histograms and definitions within tools/constants.h. Note that a few vaiables have been commented out, such as pn aand neutrino energy. This is because they were problem variables and never made it into the final analysis plots.

## Systematics
THIS FOLDER PRODUCES SYSTEMATICS FOR THE EVENT DISTRIBUTIONS ONLY. 

- make_plots.C: Produces covariance matrices for the flux and reinteraction multisims, the MicroBooNE Tune multisims, and all the MC unisims. Run using following line: root -b make_plots.C
- detVar.C: Produces covariance matrices for all of the different detector variation samples. Run using following line: root -b detVar.C
- other_detvar.C: Similar to detVar.C, but considerations for bin to bin correlations are taken into account i.e. MORE CORRECT. Run using the following:
```
root -b other_detvar.C
other_detvar s
s.main()
```
- dirt.C: Produces covariance matrices of the dirt systematic uncertainty. Dirt is treated as a unisim with one universe having 100% dirt contribution, and the second universe have 130% dirt contribution.
- statistical.C: Produces covariance matrices of the statistical uncertainty. Run using following line: root -b statistical.C
## NuWro
This folder contains code to create the event distributions using the NuWro overlay MC. The NuWro sample was generated in a similar way to the MCC9 Overlay MC, but it uses NuWro as its base MC predictor. The code produces similar plots to those created by PeLEE/analysis.C. Run the code using the following:

```
root -b nuwro_analysis.C
nuwro_analysis s
s.main()
```

## Systematics

cross section distribution found in the paper. The NuWro sample was generated in a similar way to the MCC9 Overlay MC, but it uses NuWro as its base MC predictor. It is treated exactly the
same as the overlay MC and its distribution in the plot is the denominator of the efficiency i.e. the distribbution of CC1mu2p events that adhere to the signal definition. 

