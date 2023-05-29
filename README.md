# CC2p-Analysis-Code
Ths repository contains a variety of code used to produce event distributions, MC distributions, and cross section distributions for the CC2p Analysis. 

## root_files Folder
This contains a variety of root files that were generated on the GPVMs. Further explanation is provided below

* **root_files/GCF:** CC2p signal events selected from a zombie version of GENIE that used the GCF. This model is not used in the analysis as it was confirmed that Q2 of the BNB is too low to notice these effects.
* **root_files/MEC:** These root files use different combinations of MEC event generators and if FSI are turned on or not.
* **root_files/filtered:** These are remenants from the 3 prong filter I had built. Instead of using events filtered by that module, we use the PeLEE ntuples.
* **root_files/nuwro:** In addition to the traditional Overlay MC sample, MicroBooNE also generated an overlay sample that used NuWro as the base model. These contain selected CC2p events from those samples.
* **root_files/pelee:** These are the BNB, EXT, Dirt, and Overlay CC2p events selected from the respective PeLEE ntuples. 
* **root_files/unfiltered:** I believe these were used to study raw event distributions, but I don't totally remember.

## Event Distribution Plots

### Systematics Folder
THIS FOLDER PRODUCES SYSTEMATICS FOR THE EVENT DISTRIBUTIONS ONLY. You must run the code in the following order:
 * **dirt.C:** Produces covariance matrices of the dirt systematic uncertainty. Dirt is treated as a unisim with one universe having 100% dirt contribution, and the second universe have 130% dirt contribution.
 * **statistical.C:** Produces covariance matrices of the statistical uncertainty. Run using following line: root -b statistical.C
 * **detVar.C/detVar.h:** Produces covariance matrices for all of the different detector variation samples. It treats the detector variations as unisims rather than taking the difference between each variation and the CV. Run it using the following
 ```
 # Run inside of Systematics/
 root -b detVar.C
 detVar s
 s.main()
 ```
 * **detVar_nocorrelations.C:** This differs than the above in that the difference between the CV and each variation is taken to be the systematics. This is not totally correct and the output of this code is not used in the paper. Run using following line: root -b detVar_nocorrelations.C
 * **make_plots.C:** Produces covariance matrices for the flux and reinteraction multisims, the MicroBooNE Tune multisims, and all the MC unisims. Run using following line: root -b make_plots.C

### Systematics/plotting Folder
All the code in this directory is for making plots of the various different systematics produced in the previous step. You must run it in the following order:
 * **GENIE_plot.C:** Plots the fractional uncertainty from each GENIE MC contribution and creates the total GENIE fractional uncertainty curve. Also produces the total GENIE MC correlation matrix.
 ```
 # Run inside of Systematics/plotting folder
 root -b GENIE_plot.C
 GENIE_plot s
 s.main()
 ```
 * **detvar_plot.C:** Plots the fractional uncertainty from each detector variation contribution and creates the total detector variation fractional uncertainty curve. Also produces the total detector variation correlation matrix. Run using the following:
 ```
 # Run inside of Systematics/plotting folder
 root -b detvar_plot.C
 detvar_plot s
 s.main()
 ``` 
 * **all_plot.C:** Plots the fractional uncertainty from all sources of error. root -b all_plot.C
 * **total_covariance_matrix:** Creates the total covariance matrix by adding all the covariance matrices together.
 ```
 # Run inside of Systematics/plotting folder
 root -b total_covariance_matrtix.C
 total_covariance_matrtix s
 s.main()
 ``` 
 * **total_covariance_matrix_no_stat.C:** Same as above, but doesn't include the statistical uncertainty for the purpose of creating the shape and normalization uncertainties.
 ```
 # Run inside of Systematics/plotting folder
 root -b total_covariance_matrtix_no_stat.C
 total_covariance_matrtix_no_stat s
 s.main()
 ``` 
 * **norm_and_shape.C:** Determines the normalization and shape uncertainties. Run using the following: 
 ```
 # Run inside of Systematics/plotting folder
 root -b norm_and_shape.C
 norm_and_shape s
 s.main()
 ``` 
 * **plotting.h:** Helpful plotting header file.

### PeLEE Folder
This folder contains code to create a variety of different plots using the output root files created from code found in the [Event Selection Repository]([url](https://github.com/ssfehlberg/CC2p-Event-Selection)). Produced plots include 1) efficiency distributions 2) X,Y, and Z coordinates of the reconstructed vertex and 3) the selected event distributions. There are also commented out blocks of code used to generate a variety of PFP plots and plots of cuts (such as the PID and track score value). 
```
# Run inside of PelEE Folder
root -b analysis.C
analysis s
s.main()
```
Users will then be given two different prompts. The first will indicate if you wish to look at the distributions using the nomial binning (0 = pelee) or with the optimized cross section binning (1 = pelee with xsec binning). The second prompt will ask which run you wish to process. You can find all the various histograms and definitions within tools/constants.h. Note that a few vaiables have been commented out, such as pn aand neutrino energy. This is because they were problem variables and never made it into the final analysis plots.
* **xsec_prep.C/xsec_prep.h**: Run this after you have run analysis.C. It produces the smearing matrices and the efficiency curves used in the cross section extraction. Run using the following:
```
# Run inside of PelEE Folder
root -b xsec_prep.C
xsec_prep s
s.main()
```

### NuWro Folder
This folder contains code to create the event distributions using the NuWro overlay MC. The NuWro sample was generated in a similar way to the MCC9 Overlay MC, but it uses NuWro as its base MC predictor. The code produces similar plots to those created by PeLEE/analysis.C. Run the code using the following:
```
# Run inside of NuWro Folder
root -b nuwro_analysis.C
nuwro_analysis s
s.main()
```

## Cross Sections Plots

### Xsec/RooUnfold
This folder is a clone of the [RooUnfold repository]([url](https://gitlab.cern.ch/RooUnfold/RooUnfold)). We make use of RooUnfold's implementation of D'Agostini to unfold the data from reconstructed space to true space. Note: the files in this folder are specific to my personal laptop i.e. individuals will not be able to run the following code without a system specific copy of RooUnfold. It should be as simple as redoing the make. 

### Xsec/Systematics
THIS FOLDER PRODUCES SYSTEMATICS FOR THE CROSS SECTIONS ONLY. 

This folder is very similar to the Event Distribution Systematics folder, but there are some differences. Since we have to use our unfolding procedure to determine systematics, this becomes computational intensive. All of the multisim samples were generated on the GPVM to help speed up the process. 

* **systematics.h:** The associated header file for the main script, systematics.C. All it is doing is loading in the appropriate class definitions for whatever systematic sample you wish to produce.
* **systematics.C:** The purpose of this script is to call the different systematic header files and then run them. I though I was being clever doing this, but it just turned into a big mess. Ultimately, for whatever systematic you want to run, you will have to uncomment out the relevant NAME.main() line and comment out the others. You will also have to make sure you have the correct header file included in systematics.h. To run this code you will type the following:
 ```
# Run inside of XSec/Systematics Rolder
root -b
gSystem->Load(“../RooUnfold/libRooUnfold.so”); # point to whereve that file lives on your system
.L systematics.C
systematics s
s.main()
```
* **dirt.h:** Class to produce dirt covariance matrices
* **detvar.h:** Class to produce the detector variation covariance matrices. Treated using a unisim approach.
* **mutlsims.h:** Produces the covariance matrices for the flux, reinteraction, and all MC uncertainties.

### Systematics/plotting Folder
All the code in this directory is for making plots of the various different systematics produced in the previous step. You must run it in the following order:

 * **GENIE_plot.C:** Plots the fractional uncertainty from each GENIE MC contribution and creates the total GENIE fractional uncertainty curve. Also produces the total GENIE MC correlation matrix.
 ```
 # Run inside of Xsec/Systematics/plotting folder
 root -b GENIE_plot.C
 GENIE_plot s
 s.main()
 ```
 * **detvar_plot.C:** Plots the fractional uncertainty from each detector variation contribution and creates the total detector variation fractional uncertainty curve. Also produces the total detector variation correlation matrix. Run using the following:
 ```
 # Run inside of Systematics/plotting folder
 root -b detvar_plot.C
 detvar_plot s
 s.main()
 ```
  * **all_plot.C:** Plots the fractional uncertainty from all sources of error. root -b all_plot.C
 * **total_correlation_matrix.C**: Meant to create the total correlation matrix. Run using the following:
 ```
 # Run inside of Xsec/Systematics/plotting folder
 root -b total_correlation_matrix.C
 total_correlation_matrtix s
 s.main()
 ``` 

 * **total_covariance_matrix:** Creates the total covariance matrix by adding all the covariance matrices together.
 ```
 # Run inside of Xsec/Systematics/plotting folder
 root -b total_covariance_matrtix.C
 total_covariance_matrtix s
 s.main()
 ``` 
 * **total_covariance_matrix_no_stat.C:** Same as above, but doesn't include the statistical uncertainty for the purpose of creating the shape and normalization uncertainties.
 ```
 # Run inside of Xsec/Systematics/plotting folder
 root -b total_covariance_matrtix_no_stat.C
 total_covariance_matrtix_no_stat s
 s.main()
 ``` 
 * **norm_and_shape.C:** Determines the normalization and shape uncertainties. Run using the following: 
 ```
 # Run inside of Xsec/Systematics/plotting folder
 root -b norm_and_shape.C
 norm_and_shape s
 s.main()
 ``` 
 * **plotting.h:** Helpful plotting header file.
 * **shared.h:** Seems similar to plotting.h, but seems to be used instead. Not sure why.
 * **shared_class.h** Class version of above code.

### XSec Folder
This is the main folder for creating the cross section distributions. 

* **xsec.C/xsec.h:** Produces the cross-section curves shown in the paper. If you want to run this code, do the following:
```
# Run inside of XSec Folder
root -b
gSystem->Load(“RooUnfold/libRooUnfold.so”); # point to wherever this file lives on your system.
.L xsec.C
xsec s
s.main()
```
 * **neutrino_flux.h:** Calculates the flux normalization factor and makes a pretty plot
 * **iterations.h:** Does the iteration test and computes the optimal number of unfolding iterations. 
 * **closure_test.h:** Ensures that our Smearing matrices are correct
 * **mc_model_comparison.h:** Computes the cross-section for each of the models and then puts them onto a single plot for comparison. These curves are then used in the cross-section plot with the unfolded data
 * **num_iterations.csv:** Provides the number of iterations for all the variables. This was hand made.
 * **one_iteration.csv:** An example of a single iteration that was used to determine the difference induced by using a different number of iterations.
 * **constants.h:** Variety of different parameters.  
* **CC2p_Chi2Calc.cxx:** Code from Afro to calculate the Chi2 between the data and different MC generators. It differs from a traditional chi2 in that it also uses the covariance matrix. The outputs were saved to a root file ("CC2p_xsec_results.root") to be used for further study.
* **nuwro.C:** This code was originally used in an effort to determine the number of iterations to unfold by comparing to the NuWro samples. Unfortunately, the underlying models between NuWro and the MicroBooNE Tune were so different, that this was not a helpful test. We instead compared the chi2 between adjacent numbers of iterations to determine the optimal number of unfolding computations. 
