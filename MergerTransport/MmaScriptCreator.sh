#!/bin/bash

BASEDIR="/home/dmoser/eclipse-workspace/TransportProject/MergerTransport/"
PARAMETER="XAShift"
TEMPLATEFILE="Transfer_08-09-20_NC_opt_D_basis_d"$PARAMETER".m"
TEMPLATEFILEG="Transfer_08-09-20_NC_opt_G_basis_d"$PARAMETER".m"
TEMPLATEFILEC="Transfer_08-09-20_NC_opt_C_basis_d"$PARAMETER".m"
TEMPLATEFILEH="Transfer_08-09-20_NC_opt_H_basis_d"$PARAMETER".m"
EXEDIRG='F:\"Program Files"\"Wolfram Research"\Mathematica\11.3\'
EXEDIRH='C:\"Program Files"\"Wolfram Research"\Mathematica\11.3\'
SLURMTEMPLATE="MmaSlurm_template.sh"
#bin array
lowbins=(1 99 107 123 139 155 179)
highbins=(98 106 122 138 154 178 253)
# b array
b=(-0.0015 -0.00075 0. 0.00075 0.0015)
# filename basis
for i in {0..4}
do
	bi[$i]=`echo "${b[$i]} * 10000 " | bc -l`
	bi[$i]=${bi[$i]%.*}
done
#echo ${bi[*]}

FOLDER="NC_opt_d"$PARAMETER
FILENAME="Transfer_08-09-20_"$FOLDER"_D"
FILENAMEG="Transfer_08-09-20_"$FOLDER"_G"
FILENAMEC="Transfer_08-09-20_"$FOLDER"_C"
FILENAMEH="Transfer_08-09-20_"$FOLDER"_H"

#create folder of parameter
mkdir $FOLDER
cd $FOLDER
cp ../$TEMPLATEFILE .
cp ../$TEMPLATEFILEG .
cp ../$TEMPLATEFILEC .
cp ../$TEMPLATEFILEH .
cp ../$SLURMTEMPLATE .

# now loop over different b values
for bindex in {0..4}
do
	#make subfolder for b
	mkdir "b"${bi[$bindex]}
	#create script file for this b
	SCRIPTFILE="b"${bi[$bindex]}"/"$FILENAME"_b"${bi[$bindex]}"_script.sh"
	SCRIPTFILEG="b"${bi[$bindex]}"/"$FILENAMEG"_b"${bi[$bindex]}"_script.ps1"
	SCRIPTFILEH="b"${bi[$bindex]}"/"$FILENAMEH"_b"${bi[$bindex]}"_script.ps1"
	SLURMMASTERSCRIPT="b"${bi[$bindex]}"/SlurmMaster.sh"
	echo "#!/bin/bash" > $SCRIPTFILE
	echo "#!/bin/bash" > $SLURMMASTERSCRIPT
	# loop over bin parts
	for i in {0..6}
	do
		# write into .m script files
		LOCALFILE="b"${bi[$bindex]}"/"$FILENAME"_b"${bi[$bindex]}"_P"$i".m"
		head -n 11 $TEMPLATEFILE > $LOCALFILE
		echo "b="${b[$bindex]}";" >> $LOCALFILE
		echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILE
		tail -n +14 $TEMPLATEFILE >> $LOCALFILE

		# same for G .m files
		LOCALFILEG="b"${bi[$bindex]}"/"$FILENAMEG"_b"${bi[$bindex]}"_P"$i".m"
		LOCALFILEGBACK=$FILENAMEG"_b"${bi[$bindex]}"_P"$i".m"
		head -n 11 $TEMPLATEFILEG > $LOCALFILEG
		echo "b="${b[$bindex]}";" >> $LOCALFILEG
		echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILEG
		tail -n +14 $TEMPLATEFILEG >> $LOCALFILEG
		
		# same for H .m files
		LOCALFILEH="b"${bi[$bindex]}"/"$FILENAMEH"_b"${bi[$bindex]}"_P"$i".m"
		LOCALFILEHBACK=$FILENAMEH"_b"${bi[$bindex]}"_P"$i".m"
		head -n 11 $TEMPLATEFILEH > $LOCALFILEH
		echo "b="${b[$bindex]}";" >> $LOCALFILEH
		echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILEH
		tail -n +14 $TEMPLATEFILEH >> $LOCALFILEH

		# same for C .m files
		LOCALFILEC="b"${bi[$bindex]}"/"$FILENAMEC"_b"${bi[$bindex]}"_P"$i".m"
		head -n 11 $TEMPLATEFILEC > $LOCALFILEC
		echo "b="${b[$bindex]}";" >> $LOCALFILEC
		echo "Bins = {"${lowbins[$i]}", "${highbins[$i]}"};" >> $LOCALFILEC
		tail -n +14 $TEMPLATEFILEC >> $LOCALFILEC
		
		# write script file adding each .m file for execution after each other
		echo "math -script "$LOCALFILE" >> b"${bi[$bindex]}".out" >> $SCRIPTFILE
		echo $EXEDIRG"math -script "$LOCALFILEGBACK" >> b"${bi[$bindex]}".out" >> $SCRIPTFILEG
		echo $EXEDIRH"math -script "$LOCALFILEHBACK" >> b"${bi[$bindex]}".out" >> $SCRIPTFILEH
		
		# SLURM file for Cluster
		LOCALSCRIPTFILEC="b"${bi[$bindex]}"/Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P"$i".sh"
		head -n 15 $SLURMTEMPLATE > $LOCALSCRIPTFILEC
		echo "PLACEHOLDERPATH="${FOLDER}"/b"${bi[$bindex]}"/" >> $LOCALSCRIPTFILEC
		tail -n +17 $SLURMTEMPLATE | head -n 4 >> $LOCALSCRIPTFILEC
		echo "math -script "$FILENAMEC"_b"${bi[$bindex]}"_P"$i".m" >> $LOCALSCRIPTFILEC
		echo "date" >> $LOCALSCRIPTFILEC 
		

	done

	# MASTER SCRIPT for Cluster. we write per localscript file a line, but manually because we want to include dependency stuff
	echo "jid1=\$(sbatch Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P0.sh)" >> $SLURMMASTERSCRIPT
	echo "jid1=\`echo \$jid1 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid2=\$(sbatch Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P1.sh)" >> $SLURMMASTERSCRIPT
	echo "jid2=\`echo \$jid2 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid3=\$(sbatch --dependency=afterany:\$jid1 Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P2.sh)" >> $SLURMMASTERSCRIPT
	echo "jid3=\`echo \$jid3 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid4=\$(sbatch --dependency=afterany:\$jid2 Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P3.sh)" >> $SLURMMASTERSCRIPT
	echo "jid4=\`echo \$jid4 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid5=\$(sbatch --dependency=afterany:\$jid3 Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P4.sh)" >> $SLURMMASTERSCRIPT
	echo "jid5=\`echo \$jid5 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid6=\$(sbatch --dependency=afterany:\$jid4 Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P5.sh)" >> $SLURMMASTERSCRIPT
	echo "jid6=\`echo \$jid6 | cut -d' ' -f 4\`" >> $SLURMMASTERSCRIPT
	echo "jid7=\$(sbatch --dependency=afterany:\$jid5 Slurm_"$FILENAMEC"_b"${bi[$bindex]}"_P6.sh)" >> $SLURMMASTERSCRIPT
	
	
done
