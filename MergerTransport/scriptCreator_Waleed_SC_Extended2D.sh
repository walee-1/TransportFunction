#!/bin/bash

if [[ $# -eq 0 ]]
then
	echo "Please provide Parameter name"
	exit 1
fi

if [[ $# -ge 1 ]]
then
	PARAMETER=$1
	echo "Parameter set to $PARAMETER"
fi

#offset
offset=( -0.00015 0.00005 0.00015 )

BASEDIR="/home/waleed/Documents/Wolfram_Mathematica/TransportFunction/MergerTransport"

#Paths are cluster/C, work/D, Gertrud/G
PATHC="/users/waleed.khalid/Mma/"
PATHH="/home/dmoser/Waleed/"
PATHD="/home/waleed/Documents/Wolfram_Mathematica/TransportFunction"
PATHG="C:/Users/smi/Desktop/WaleedTransport/"
TEMPLATEFILE="Transfer_08-09-20_SC_opt_basis_d"$PARAMETER"_2D.m"
EXEDIRG='F:\"Program Files"\"Wolfram Research"\Mathematica\11.3\'

SLURMTEMPLATE="MmaSlurm_template.sh"

if [[ ! -f $TEMPLATEFILE ]]
then
	echo "File does not exist. Check again"
	exit 1
fi

KERNELSD=4
KERNELSG=6
KERNELSC=8
KERNELSH=4
#bin array
lowbins=(1 50 99 107 123 139 155 179)
highbins=(49 98 106 122 138 154 178 264)
# a array
a=( -0.106 -0.1055 -0.105 -0.1045 -0.104 )
# filename basis
for i in ${!a[@]}
do
	ai[$i]=`echo "${a[$i]} * 10000 " | bc -l`
	ai[$i]=${ai[$i]%.*}
done

for i in ${!offset[@]}
do
    offseti[$i]=`echo "${offset[$i]} * 100000 " | bc -l`
	offseti[$i]=${offseti[$i]%.*}
done

#echo ${bi[*]}

FOLDER="SC_Prot_opt_d"$PARAMETER"_2D"
FILENAME="Transfer_08-09-20_"$FOLDER"_D"
FILENAMEG="Transfer_08-09-20_"$FOLDER"_G"
FILENAMEC="Transfer_08-09-20_"$FOLDER"_C"
FILENAMEH="Transfer_08-09-20_"$FOLDER"_H"

#create folder of parameter
mkdir -p $FOLDER
cd $FOLDER
cp ../$TEMPLATEFILE .
#cp ../$TEMPLATEFILEG .
#cp ../$TEMPLATEFILEC .
#cp ../$TEMPLATEFILEH .
cp ../$SLURMTEMPLATE .

writer_Function(){
    option=$1
    if [[ $option == "D" ]] 
    then
        LOCALFILE=$createdFolder"/"$FILENAME"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m"
        TEMPPATH=$PATHD
        KERNELS=$KERNELSD
		LOCALFILEDBACK=$FILENAME"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m"
		echo "math -script "$LOCALFILEDBACK" >> a"${ai[$bindex]}"_"${offseti[$bindex2]}".out" >> $SCRIPTFILE
	elif [[ $option == "H" ]] 
    then
        LOCALFILE=$createdFolder"/"$FILENAMEH"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m"
        TEMPPATH=$PATHH
        KERNELS=$KERNELSH
		LOCALFILEHBACK=$FILENAMEH"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m"
		echo "math -script "$LOCALFILEHBACK" >> a"${ai[$bindex]}"_"${offseti[$bindex2]}".out" >> $SCRIPTFILEH
    elif [[ $option == "G" ]] 
    then
        LOCALFILE=$createdFolder"/"$FILENAMEG"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m"
        TEMPPATH=$PATHG
        KERNELS=$KERNELSG
        LOCALFILEGBACK=$FILENAMEG"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m"
		echo $EXEDIRG"math -script "$LOCALFILEGBACK" >> a"${ai[$bindex]}"_"${offseti[$bindex2]}".out" >> $SCRIPTFILEG
    elif [[ $option == "C" ]] 
    then
        LOCALFILE=$createdFolder"/"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m"
        TEMPPATH=$PATHC
        KERNELS=$KERNELSC
		 # SLURM file for Cluster
	    LOCALSCRIPTFILEC="a"${ai[$bindex]}"_"${offseti[$bindex2]}"/Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".sh"
	    head -n 15 $SLURMTEMPLATE > $LOCALSCRIPTFILEC
	    echo "PLACEHOLDERPATH="${FOLDER}"/a"${ai[$bindex]}"_"${offseti[$bindex2]}"/" >> $LOCALSCRIPTFILEC
	    tail -n +17 $SLURMTEMPLATE | head -n 4 >> $LOCALSCRIPTFILEC
	    echo "math -script "$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P"$i".m" >> $LOCALSCRIPTFILEC
	    echo "date" >> $LOCALSCRIPTFILEC 
    fi
    cat $TEMPLATEFILE>$LOCALFILE
    sed -i "s+^kernels=.*$+kernels=$KERNELS+" $LOCALFILE
    sed -i "s+^SetDirectory.*$+SetDirectory[\"$TEMPPATH\"];+" $LOCALFILE
    sed -i "s+^a=.*$+a=${a[$bindex]};+" $LOCALFILE
    sed -i "s+^offset=.*$+offset=${offset[$bindex2]};+" $LOCALFILE
    sed -i "s+^Bins .*$+Bins = {${lowbins[$i]}, ${highbins[$i]}};+" $LOCALFILE
    sed -i "s+^Export.*$+Export[\"MergerTransport/$FOLDER/a${ai[$bindex]}_${offseti[$bindex2]}/TransferResult_08-09-20_SC_opt_$option\_d$PARAMETER\_a${ai[$bindex]}_${offseti[$bindex2]}_\"<>ToString[Bins[[1]]]<>\"-\"<>ToString[Bins[[2]]]<>\".txt\",BinwYShiftPrec44OriginalBinsb0NewBGrad,\"Table\"];+" $LOCALFILE

}

# now loop over different b values
for bindex in ${!a[@]}
do
    for bindex2  in ${!offset[@]}
    do

        #make subfolder for b
        createdFolder="a"${ai[$bindex]}"_"${offseti[$bindex2]}
        mkdir -p $createdFolder
        
        #create script file for this b
        SCRIPTFILE=$createdFolder"/"$FILENAME"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_script.sh"
        SCRIPTFILEG=$createdFolder"/"$FILENAMEG"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_script.ps1"
        if [ -f $SCRIPTFILEG ]; then
            rm $SCRIPTFILEG
        fi
        SCRIPTFILEH=$createdFolder"/"$FILENAMEH"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_script.sh"
        SLURMMASTERSCRIPT=$createdFolder"/SlurmMaster.sh"
        echo "#!/bin/bash" > $SCRIPTFILE
        echo "#!/bin/bash" > $SCRIPTFILEH
        echo "#!/bin/bash" > $SLURMMASTERSCRIPT
        # loop over bin parts
        for i in ${!lowbins[@]}
        do
            writer_Function "D"
            writer_Function "G"
            writer_Function "C"
            writer_Function "H"
        done
	# MASTER SCRIPT for Cluster. we write per localscript file a line, but manually because we want to include dependency stuff
	echo "jid1=\$(sbatch Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P0.sh)" >> $SLURMMASTERSCRIPT
	echo "jid1=\`echo \$jid1 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid2=\$(sbatch Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P1.sh)" >> $SLURMMASTERSCRIPT
	echo "jid2=\`echo \$jid2 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid3=\$(sbatch --dependency=afterany:\$jid1 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P2.sh)" >> $SLURMMASTERSCRIPT
	echo "jid3=\`echo \$jid3 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid4=\$(sbatch --dependency=afterany:\$jid2 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P3.sh)" >> $SLURMMASTERSCRIPT
	echo "jid4=\`echo \$jid4 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid5=\$(sbatch --dependency=afterany:\$jid3 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P4.sh)" >> $SLURMMASTERSCRIPT
	echo "jid5=\`echo \$jid5 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid6=\$(sbatch --dependency=afterany:\$jid4 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P5.sh)" >> $SLURMMASTERSCRIPT
	echo "jid6=\`echo \$jid6 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid7=\$(sbatch --dependency=afterany:\$jid5 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P6.sh)" >> $SLURMMASTERSCRIPT
	echo "jid7=\`echo \$jid7 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid8=\$(sbatch --dependency=afterany:\$jid6 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_"${offseti[$bindex2]}"_P7.sh)" >> $SLURMMASTERSCRIPT
	
	done
done

echo "Process Completed. Please check folder: $FOLDER"
echo "Cluster Copy Code:"
echo "scp -r $FOLDER waleed.khalid@cbe.vbc.ac.at:/users/waleed.khalid/Mma/MergerTransport/"