#!/bin/bash

BASEDIR="/home/waleed/Documents/Wolfram_Mathematica/TransportFunction/MergerTransport"
PARAMETER="a0"

#Paths are cluster/C, work/D, Gertrud/G
PATHC="/users/waleed.khalid/Mma/"
PATHD="/home/waleed/Documents/Wolfram_Mathematica/TransportFunction"
PATHG="C:/Users/smi/Desktop/WaleedTransport/"
TEMPLATEFILE="Transfer_08-09-20_SC_opt_basis_d"$PARAMETER".m"
#TEMPLATEFILEG="Transfer_08-09-20_SC_opt_G_basis_d"$PARAMETER".m"
#TEMPLATEFILEC="Transfer_08-09-20_SC_opt_C_basis_d"$PARAMETER".m"
#TEMPLATEFILEH="Transfer_08-09-20_NC_opt_H_basis_d"$PARAMETER".m"
#EXEDIRS=('F:\"Program Files"\"Wolfram Research"\Mathematica\11.3\')
EXEDIRG='F:\"Program Files"\"Wolfram Research"\Mathematica\11.3\'
#EXEDIRH='C:\"Program Files"\"Wolfram Research"\Mathematica\11.3\'
SLURMTEMPLATE="MmaSlurm_template.sh"
#bin array
lowbins=(1 99 107 123 139 155 179)
highbins=(98 106 122 138 154 178 264)
# a array
a=(-0.106 -0.1055 -0.105 -0.1045 -0.104 )
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
#FILENAMEH="Transfer_08-09-20_"$FOLDER"_H"

#create folder of parameter
mkdir -p $FOLDER
cd $FOLDER
cp ../$TEMPLATEFILE .
cp ../$TEMPLATEFILEG .
cp ../$TEMPLATEFILEC .
#cp ../$TEMPLATEFILEH .
cp ../$SLURMTEMPLATE .

writer_Function(){
    option=$1
    if [ $option -eq "D" ] then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAME"_a"${ai[$bindex]}"_P"$i".m"
    elif [ $option -eq "G" ] then
        LOCALFILE="a"${ai[$bindex]}"/"$FILENAMEG"_a"${ai[$bindex]}"_P"$i".m"
    elif [ $option -eq "C" ] then


    fi
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
		# write into .m script files
		LOCALFILE="a"${ai[$bindex]}"/"$FILENAME"_a"${ai[$bindex]}"_P"$i".m"
        cat $TEMPLATEFILE>$LOCALFILE
        sed -i "s+^SetDirectory.*$+SetDirectory[\"$PATHD\"];+" $LOCALFILE
        sed -i "s+^a=.*$+a=${a[$bindex]};+" $LOCALFILE
        sed -i "s+^Bins .*$+Bins = {${lowbins[$i]}, ${highbins[$i]}};+" $LOCALFILE
        sed -i "s+^Export.$+\Export[\"MergerTransport/SC_Prot_opt_d"$PARAMETER"\"<>ToString[IntegerPart[a*10000]]<>\"/TransferResult_SC_Prot_opt_C_d"$PARAMETER"\"<>ToString[IntegerPart[a*10000]]<>\"]"
		#head -n 11 $TEMPLATEFILE > $LOCALFILE
		#echo "a="${a[$bindex]}";" >> $LOCALFILE
		#echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILE
		#tail -n +14 $TEMPLATEFILE >> $LOCALFILE

		# same for G .m files
		LOCALFILEG="a"${ai[$bindex]}"/"$FILENAMEG"_a"${ai[$bindex]}"_P"$i".m"
		LOCALFILEGBACK=$FILENAMEG"_a"${ai[$bindex]}"_P"$i".m"
		cat $TEMPLATEFILEG > $LOCALFILEG
        sed -i "s+^SetDirectory.*$+SetDirectory[\"$PATHG\"];+" $LOCALFILEG
        sed -i "s+^a=.*$+a=${a[$bindex]};+" $LOCALFILEG
        sed -i "s+^Bins .*$+Bins = {${lowbins[$i]}, ${highbins[$i]}};+" $LOCALFILEG
		#echo "a="${a[$bindex]}";" >> $LOCALFILEG
		#echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILEG

		
		# same for H .m files
#		LOCALFILEH="b"${bi[$bindex]}"/"$FILENAMEH"_b"${bi[$bindex]}"_P"$i".m"
#		LOCALFILEHBACK=$FILENAMEH"_b"${bi[$bindex]}"_P"$i".m"
#		head -n 11 $TEMPLATEFILEH > $LOCALFILEH
#		echo "b="${b[$bindex]}";" >> $LOCALFILEH
#		echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILEH
#		tail -n +14 $TEMPLATEFILEH >> $LOCALFILEH

		# same for C .m files
		LOCALFILEC="a"${ai[$bindex]}"/"$FILENAMEC"_a"${ai[$bindex]}"_P"$i".m"
		cat $TEMPLATEFILEC > $LOCALFILEC
        sed -i "s+^SetDirectory.*$+SetDirectory[\"$PATHC\"];+" $LOCALFILEC
        sed -i "s+^a=.*$+a=${a[$bindex]};+" $LOCALFILEC
        sed -i "s+^Bins .*$+Bins = {${lowbins[$i]}, ${highbins[$i]}};+" $LOCALFILEC
		#echo "a="${a[$aindex]}";" >> $LOCALFILEC
		#echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILEC
		#tail -n +14 $TEMPLATEFILEC >> $LOCALFILEC
		
		# write script file adding each .m file for execution after each other
		echo "math -script "$LOCALFILE" >> a"${ai[$bindex]}".out" >> $SCRIPTFILE
		echo $EXEDIRG"math -script "$LOCALFILEGBACK" >> a"${ai[$bindex]}".out" >> $SCRIPTFILEG
#		echo $EXEDIRH"math -script "$LOCALFILEHBACK" >> a"${ai[$bindex]}".out" >> $SCRIPTFILEH
		
		# SLURM file for Cluster
		LOCALSCRIPTFILEC="a"${ai[$bindex]}"/Slurm_"$FILENAMEC"_a"${ai[$bindex]}"_P"$i".sh"
		head -n 15 $SLURMTEMPLATE > $LOCALSCRIPTFILEC
		echo "PLACEHOLDERPATH="${FOLDER}"/a"${ai[$bindex]}"/" >> $LOCALSCRIPTFILEC
		tail -n +17 $SLURMTEMPLATE | head -n 4 >> $LOCALSCRIPTFILEC
		echo "math -script "$FILENAMEC"_a"${ai[$bindex]}"_P"$i".m" >> $LOCALSCRIPTFILEC
		echo "date" >> $LOCALSCRIPTFILEC 
		

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
