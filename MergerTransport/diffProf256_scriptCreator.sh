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

#Paths are cluster/C, work/D, Gertrud/G
PATHC="/users/waleed.khalid/Mma/"
PATHH="/home/dmoser/Waleed/"
PATHD="/home/waleed/Documents/Wolfram_Mathematica/TransportFunction"
PATHG="C:/Users/smi/Desktop/WaleedTransport/" 
TEMPLATEFILE="Transfer_08-09-20_SC_diffProf_d"$PARAMETER".m"
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
lowbins=( 1 34 67 81 100 111 121 136 151 161 171 181 191 201 210 226 232 )
highbins=( 33 66 80 99 110 120 135 150 160 170 180 190 200 209 225 231 256 )
# a array
a=( 0.105 -0.106 -0.1055 -0.105 -0.1045 -0.104 )
# filename basis
for i in ${!a[@]}
do
	ai[$i]=`echo "${a[$i]} * 10000 " | bc -l`
	ai[$i]=${ai[$i]%.*}
done
#echo ${bi[*]}

FOLDER="SC_Prot_diffProf_d"$PARAMETER
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
		LOCALFILEDBACK=$FILENAME"_a"${ai[$bindex]}"_P"$i".m"
		echo "math -script "$LOCALFILEDBACK" >> a"${ai[$bindex]}".out" >> $SCRIPTFILE
		echo "echo Finished: $LOCALFILEDBACK $(date)" >> $SCRIPTFILE
	elif [[ $option == "H" ]] 
    then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAMEH"_a"${ai[$bindex]}"_P"$i".m"
        TEMPPATH=$PATHH
        KERNELS=$KERNELSH
		LOCALFILEHBACK=$FILENAMEH"_a"${ai[$bindex]}"_P"$i".m"
		echo "math -script "$LOCALFILEHBACK" >> a"${ai[$bindex]}".out" >> $SCRIPTFILEH
		echo "echo Finished: $LOCALFILEHBACK $(date)" >> $SCRIPTFILEH
    elif [[ $option == "G" ]] 
    then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAMEG"_a"${ai[$bindex]}"_P"$i".m"
        TEMPPATH=$PATHG
        KERNELS=$KERNELSG
        LOCALFILEGBACK=$FILENAMEG"_a"${ai[$bindex]}"_P"$i".m"
		echo $EXEDIRG"math -script "$LOCALFILEGBACK" >> a"${ai[$bindex]}".out" >> $SCRIPTFILEG
		echo "echo Finished: $LOCALFILEGBACK $(date)" >> $SCRIPTFILEG
    elif [[ $option == "C" ]] 
    then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAMEC"_a"${ai[$bindex]}"_P"$i".m"
        TEMPPATH=$PATHC
        KERNELS=$KERNELSC
		 # SLURM file for Cluster
	    LOCALSCRIPTFILEC="a"${ai[$bindex]}"/Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P"$i".sh"
	    head -n 15 $SLURMTEMPLATE > $LOCALSCRIPTFILEC
		sed -i "s+#SBATCH --job-name=Mma_TransP1+#SBATCH --job-name=${PARAMETER}_P$i+" $LOCALSCRIPTFILEC
	    echo "PLACEHOLDERPATH="${FOLDER}"/a"${ai[$bindex]}"/" >> $LOCALSCRIPTFILEC
	    tail -n +17 $SLURMTEMPLATE | head -n 4 >> $LOCALSCRIPTFILEC
	    echo "math -script "$FILENAMEC"_a"${ai[$bindex]}"_P"$i".m" >> $LOCALSCRIPTFILEC
	    echo "date" >> $LOCALSCRIPTFILEC 
    fi
    cat $TEMPLATEFILE>$LOCALFILE
    sed -i "s+^kernels=.*$+kernels=$KERNELS;+" $LOCALFILE
    sed -i "s+^SetDirectory.*$+SetDirectory[\"$TEMPPATH\"];+" $LOCALFILE
    sed -i "s+^a=.*$+a=${a[$bindex]};+" $LOCALFILE
    sed -i "s+^Bins .*$+Bins = {${lowbins[$i]}, ${highbins[$i]}};+" $LOCALFILE
    sed -i "s+^Export.*$+Export[\"MergerTransport/SC_Prot_diffProf_d$PARAMETER/a\"<>ToString[IntegerPart[a*10000]]<>\"/TransferResult_08-09-20_SC_diffProf_$option\_d$PARAMETER\_a\"<>ToString[IntegerPart[a*10000]]<>\"_\"<>ToString[Bins[[1]]]<>\"-\"<>ToString[Bins[[2]]]<>\".txt\",BinwYShiftPrec44OriginalBinsb0NewBGrad,\"Table\"];+" $LOCALFILE

}

# now loop over different b values
for bindex in ${!a[@]}
do
	#make subfolder for b
	mkdir -p "a"${ai[$bindex]}
	#create script file for this b
	SCRIPTFILE="a"${ai[$bindex]}"/"$FILENAME"_a"${ai[$bindex]}"_script.sh"
	SCRIPTFILEG="a"${ai[$bindex]}"/"$FILENAMEG"_a"${ai[$bindex]}"_script.ps1"
	if [ -f $SCRIPTFILEG ]; then
		rm $SCRIPTFILEG
	fi
	SCRIPTFILEH="a"${ai[$bindex]}"/"$FILENAMEH"_a"${ai[$bindex]}"_script.sh"
	SLURMMASTERSCRIPT="a"${ai[$bindex]}"/SlurmMaster.sh"
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


	binNum=${#lowbins[@]}
	remJIDs="$((binNum-2))"

	# MASTER SCRIPT for Cluster. we write per localscript file a line, but manually because we want to include dependency stuff
	for(( i=1;i<=$binNum;i++ ));
	do
        pNo="$((i-1))"
        prevJID="$((i-2))"
        if [[ $i -le 2 ]]
        then
		    echo "jid"${i}"=\$(sbatch Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P"${pNo}".sh)" >> $SLURMMASTERSCRIPT
		    echo "jid"${i}"=\`echo \$jid"${i}" | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
        elif [[ $i -lt $binNum ]]
        then
            echo "jid"${i}"=\$(sbatch --dependency=afterany:\$jid"${prevJID}" Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P"${pNo}".sh)" >> $SLURMMASTERSCRIPT
	        echo "jid"${i}"=\`echo \$jid"${i}" | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
        else
            echo "jid"${i}"=\$(sbatch --dependency=afterany:\$jid"${prevJID}" Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P"${pNo}".sh)" >> $SLURMMASTERSCRIPT
        fi
            
	done


# 	# MASTER SCRIPT for Cluster. we write per localscript file a line, but manually because we want to include dependency stuff
# 	echo "jid1=\$(sbatch Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P0.sh)" >> $SLURMMASTERSCRIPT
# 	echo "jid1=\`echo \$jid1 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
# 	echo "jid2=\$(sbatch Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P1.sh)" >> $SLURMMASTERSCRIPT
# 	echo "jid2=\`echo \$jid2 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
# 	echo "jid3=\$(sbatch --dependency=afterany:\$jid1 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P2.sh)" >> $SLURMMASTERSCRIPT
# 	echo "jid3=\`echo \$jid3 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
# 	echo "jid4=\$(sbatch --dependency=afterany:\$jid2 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P3.sh)" >> $SLURMMASTERSCRIPT
# 	echo "jid4=\`echo \$jid4 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
# 	echo "jid5=\$(sbatch --dependency=afterany:\$jid3 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P4.sh)" >> $SLURMMASTERSCRIPT
# 	echo "jid5=\`echo \$jid5 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
# 	echo "jid6=\$(sbatch --dependency=afterany:\$jid4 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P5.sh)" >> $SLURMMASTERSCRIPT
# 	echo "jid6=\`echo \$jid6 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
# 	echo "jid7=\$(sbatch --dependency=afterany:\$jid5 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P6.sh)" >> $SLURMMASTERSCRIPT
# 	echo "jid7=\`echo \$jid7 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
# 	echo "jid8=\$(sbatch --dependency=afterany:\$jid6 Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P7.sh)" >> $SLURMMASTERSCRIPT
	
	
done

echo "Process Completed. Please check folder: $FOLDER"
echo "Cluster Copy Code:"
echo "scp -r $FOLDER waleed.khalid@cbe.vbc.ac.at:/users/waleed.khalid/Mma/MergerTransport/"