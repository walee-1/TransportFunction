#!/bin/bash


binNum=8
if [[ $# -lt 2 ]]
then
	echo "Usage: DirName a-array"
	exit 1
fi

if [[ $# -ge 2 ]]
then
	dirName=$1; shift
    a=( "$@" )
	echo "Parameter set to $dirName"
    echo "sets of bins set to $binNum"
fi
cd $dirName

FILENAMEC="Transfer_08-09-20_"$dirName"_C"

scriptFile="masterScript.sh"

echo "#!/bin/bash" > $scriptFile
echo "#SBATCH --nodes=1" >>$scriptFile
echo "#SBATCH --ntasks=1" >>$scriptFile
echo "#SBATCH --time=00:01:00" >>$scriptFile
echo "#SBATCH --output=output.txt">>$scriptFile
j=1
for bindex in ${!a[@]}
do
        
        if [[ $bindex -eq 0 ]];then
            echo "cd a${a[$bindex]}" >>$scriptFile
        else
            echo "cd ../a${a[$bindex]}" >>$scriptFile
        fi
        
        for(( i=1;i<=$binNum;i++ ));
        do
            pNo="$((i-1))"
            prevJID="$((j-2))"
            if [[ $j -le 2 ]]
            then
                echo "jid"${j}"=\$(sbatch Slurm_"$FILENAMEC"_a"${a[$bindex]}"_P"${pNo}".sh)" >> $scriptFile
                echo "jid"${j}"=\`echo \$jid"${j}" | cut -d' ' -f 4\`" >> $scriptFile
            else    
                echo "jid"${j}"=\$(sbatch --dependency=afterany:\$jid"${prevJID}" Slurm_"$FILENAMEC"_a"${a[$bindex]}"_P"${pNo}".sh)" >> $scriptFile
                echo "jid"${j}"=\`echo \$jid"${j}" | cut -d' ' -f 4\`" >> $scriptFile
            fi
            j="$((j+1))"
        done

done
cd ../
./PCServer_CopyScript.sh $dirName
