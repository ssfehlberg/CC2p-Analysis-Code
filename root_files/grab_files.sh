#!/bin/bash

echo "Grabbing Run 1 Files"
cd root_files/pelee/Run1/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/PeLEE_ntuples/root_files/Run1/*.root .
cd -
cd root_files/nuwro/Run1/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/NuWro/Run1/*.root .
cd -

echo "Grabbing Run 2 Files"
cd root_files/pelee/Run2/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/PeLEE_ntuples/root_files/Run2/*.root .
cd -
cd root_files/nuwro/Run2/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/NuWro/Run2/*.root .
cd -

echo "Grabbing Run 3 Files"
cd root_files/pelee/Run3/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/PeLEE_ntuples/root_files/Run3/*.root .
cd -
cd root_files/nuwro/Run3/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/NuWro/Run3/*.root .
cd -

echo "Grabbing Run All Files"
cd root_files/pelee/Run_all/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/PeLEE_ntuples/root_files/Run_all/*.root .
cd -
cd root_files/nuwro/Run_all/
scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/2proton_Pandora/NuWro/Run_all/*.root .
cd -
