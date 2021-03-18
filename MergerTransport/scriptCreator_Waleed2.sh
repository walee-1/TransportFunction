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


BASEDIR="/home/waleed/Documents/Wolfram_Mathematica/TransportFunction/MergerTransport"
PARAMETER="a0"

#Paths are cluster/C, work/D, Gertrud/G
PATHC="/users/waleed.khalid/Mma/"
PATHH="/home/dmoser/Waleed/"
PATHD="/home/waleed/Documents/Wolfram_Mathematica/TransportFunction"
PATHG="C:/Users/smi/Desktop/WaleedTransport/"
TEMPLATEFILE="Transfer_08-09-20_SC_opt_basis_d"$PARAMETER".m"
EXEDIRG='F:\"Program Files"\"Wolfram Research"\Mathematica\11.3\'

SLURMTEMPLATE="MmaSlurm_template.sh"

if [[ ! -f $TEMPLATEFILE ]]
then
	echo "File does not exist. Check again"
	exit 1
fi

KERNELSD=6
KERNELSG=6
KERNELSC=8
KERNELSH=4
#bin array
lowbins=(1 99 107 123 139 155 179)
highbins=(98 106 122 138 154 178 264)
# a array
a=( -0.106 -0.1055 -0.105 -0.1045 -0.104 )
# filename basis
for i in {0..4}
do
	ai[$i]=`echo "${a[$i]} * 10000 " | bc -l`
	ai[$i]=${ai[$i]%.*}
done
#echo ${bi[*]}

FOLDER="SC_Prot_opt_d"$PARAMETER
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
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAME"_a"${ai[$bindex]}"_P"$i".m"
        TEMPPATH=$PATHD
        KERNELS=$KERNELSD
		echo "math -script "$LOCALFILE" >> a"${ai[$bindex]}".out" >> $SCRIPTFILE
	elif [[ $option == "H" ]] 
    then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAMEH"_a"${ai[$bindex]}"_P"$i".m"
        TEMPPATH=$PATHH
        KERNELS=$KERNELSD
		echo "math -script "$LOCALFILE" >> a"${ai[$bindex]}".out" >> $SCRIPTFILE
    elif [[ $option == "G" ]] 
    then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAMEG"_a"${ai[$bindex]}"_P"$i".m"
        TEMPPATH=$PATHG
        KERNELS=$KERNELSG
        LOCALFILEGBACK=$FILENAMEG"_a"${ai[$bindex]}"_P"$i".m"
		echo $EXEDIRG"math -script "$LOCALFILEGBACK" >> a"${ai[$bindex]}".out" >> $SCRIPTFILEG
    elif [[ $option == "C" ]] 
    then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAMEC"_a"${ai[$bindex]}"_P"$i".m"
        TEMPPATH=$PATHC
        KERNELS=$KERNELSC
		 # SLURM file for Cluster
	    LOCALSCRIPTFILEC="a"${ai[$bindex]}"/Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P"$i".sh"
	    head -n 15 $SLURMTEMPLATE > $LOCALSCRIPTFILEC
	    echo "PLACEHOLDERPATH="${FOLDER}"/a"${ai[$bindex]}"/" >> $LOCALSCRIPTFILEC
	    tail -n +17 $SLURMTEMPLATE | head -n 4 >> $LOCALSCRIPTFILEC
	    echo "math -script "$FILENAMEC"_a"${ai[$bindex]}"_P"$i".m" >> $LOCALSCRIPTFILEC
	    echo "date" >> $LOCALSCRIPTFILEC 
    fi
    cat $TEMPLATEFILE>$LOCALFILE
    sed -i "s+^kernels=.*$+kernels=$KERNELS+" $LOCALFILE
    sed -i "s+^SetDirectory.*$+SetDirectory[\"$TEMPPATH\"];+" $LOCALFILE
    sed -i "s+^a=.*$+a=${a[$bindex]};+" $LOCALFILE
    sed -i "s+^Bins .*$+Bins = {${lowbins[$i]}, ${highbins[$i]}};+" $LOCALFILE
    sed -i "s+^Export.*$+Export[\"MergerTransport/SC_Prot_opt_d$PARAMETER/a\"<>ToString[IntegerPart[a*10000]]<>\"/TransferResult_08-09-20_SC_opt_$option\_d$PARAMETER\_a\"<>ToString[IntegerPart[a*10000]]<>\"_\"<>ToString[Bins[[1]]]<>\"-\"<>ToString[Bins[[2]]]<>\".txt\",BinwYShiftPrec44OriginalBinsb0NewBGrad,\"Table\"];+" $LOCALFILE

}

# now loop over different b values
for bindex in {0..4}
do
	#make subfolder for b
	mkdir -p "a"${ai[$bindex]}
	#create script file for this b
	SCRIPTFILE="a"${ai[$bindex]}"/"$FILENAME"_a"${ai[$bindex]}"_script.sh"
	SCRIPTFILEG="a"${ai[$bindex]}"/"$FILENAMEG"_a"${ai[$bindex]}"_script.ps1"
#	SCRIPTFILEH="a"${ai[$bindex]}"/"$FILENAMEH"_a"${ai[$bindex]}"_script.ps1"
	SLURMMASTERSCRIPT="a"${ai[$bindex]}"/SlurmMaster.sh"
	echo "#!/bin/bash" > $SCRIPTFILE
	echo "#!/bin/bash" > $SLURMMASTERSCRIPT
	# loop over bin parts
	for i in {0..6}
	do
        writer_Function "D"
        writer_Function "G"
        writer_Function "C"
		writer_Function "H"
	done

	# MASTER SCRIPT for Cluster. we write per localscript file a line, but manually because we want to include dependency stuff
	echo "jid1=\$(sbatch Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P0.sh)" >> $SLURMMASTERSCRIPT
	echo "jid1=\`echo \$jid1 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid2=\$(sbatch Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P1.sh)" >> $SLURMMASTERSCRIPT
	echo "jid2=\`echo \$jid2 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid3=\$(sbatch --dependency=afterany:\$jid1 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P2.sh)" >> $SLURMMASTERSCRIPT
	echo "jid3=\`echo \$jid3 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid4=\$(sbatch --dependency=afterany:\$jid2 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P3.sh)" >> $SLURMMASTERSCRIPT
	echo "jid4=\`echo \$jid4 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid5=\$(sbatch --dependency=afterany:\$jid3 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P4.sh)" >> $SLURMMASTERSCRIPT
	echo "jid5=\`echo \$jid5 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid6=\$(sbatch --dependency=afterany:\$jid4 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P5.sh)" >> $SLURMMASTERSCRIPT
	echo "jid6=\`echo \$jid6 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid7=\$(sbatch --dependency=afterany:\$jid5 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P6.sh)" >> $SLURMMASTERSCRIPT
	
	
done

echo "Process Completed. Please check folder: $FOLDER"
echo "Cluster Copy Code:"
echo "scp -r $FOLDER waleed.khalid@cbe.vbc.ac.at:/users/waleed.khalid/Mma/MergerTransport/"