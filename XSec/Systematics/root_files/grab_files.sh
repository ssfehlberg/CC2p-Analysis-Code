#!/bin/bash
input="directories.list"
while IFS= read -r line
do

    #if [ "$line" == "detVar" ]; then
#	cd "$line"
#	echo "Inside of $line"
#	echo "Grabbing /uboone/data/users/sfehlber/Systematics/$line/*.root"
#	scp sfehlber@uboonegpvm03.fnal.gov:/uboone/data/users/sfehlber/Systematics/"$line"/*.root .
#	cd -
 #   else
	cd "$line"
	echo "Inside of $line"
	echo "Grabbing /uboone/data/users/sfehlber/Systematics/$line/histograms_pelee_xsec.root"
#	scp sfehlber@uboonegpvm03.fnal.gov:/uboone/data/users/sfehlber/Systematics/"$line"/*.root .
	scp sfehlber@uboonegpvm03.fnal.gov:/uboone/app/users/sfehlber/Analysis/CC2p/Analysis/systematics/root_files/"$line"/systematics.root .
	cd -	#  cd "$line"
 #   fi

done < "$input"
      
#code to make the directories
#while IFS= read -r line
#do
#    mkdir "$line"
#done < "$input"
