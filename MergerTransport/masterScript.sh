#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00
#SBATCH --output=output.txt
cd a-1040
jid1=$(sbatch Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P0.sh)
jid1=`echo $jid1 | cut -d' ' -f 4`
jid2=$(sbatch Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P1.sh)
jid2=`echo $jid2 | cut -d' ' -f 4`
jid3=$(sbatch --dependency=afterany:$jid1 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P2.sh)
jid3=`echo $jid3 | cut -d' ' -f 4`
jid4=$(sbatch --dependency=afterany:$jid2 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P3.sh)
jid4=`echo $jid4 | cut -d' ' -f 4`
jid5=$(sbatch --dependency=afterany:$jid3 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P4.sh)
jid5=`echo $jid5 | cut -d' ' -f 4`
jid6=$(sbatch --dependency=afterany:$jid4 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P5.sh)
jid6=`echo $jid6 | cut -d' ' -f 4`
jid7=$(sbatch --dependency=afterany:$jid5 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P6.sh)
jid7=`echo $jid7 | cut -d' ' -f 4`
jid8=$(sbatch --dependency=afterany:$jid6 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1040_P7.sh)
jid8=`echo $jid8 | cut -d' ' -f 4`
cd ../a-1045
jid9=$(sbatch --dependency=afterany:$jid7 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P0.sh)
jid9=`echo $jid9 | cut -d' ' -f 4`
jid10=$(sbatch --dependency=afterany:$jid8 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P1.sh)
jid10=`echo $jid10 | cut -d' ' -f 4`
jid11=$(sbatch --dependency=afterany:$jid9 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P2.sh)
jid11=`echo $jid11 | cut -d' ' -f 4`
jid12=$(sbatch --dependency=afterany:$jid10 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P3.sh)
jid12=`echo $jid12 | cut -d' ' -f 4`
jid13=$(sbatch --dependency=afterany:$jid11 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P4.sh)
jid13=`echo $jid13 | cut -d' ' -f 4`
jid14=$(sbatch --dependency=afterany:$jid12 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P5.sh)
jid14=`echo $jid14 | cut -d' ' -f 4`
jid15=$(sbatch --dependency=afterany:$jid13 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P6.sh)
jid15=`echo $jid15 | cut -d' ' -f 4`
jid16=$(sbatch --dependency=afterany:$jid14 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1045_P7.sh)
jid16=`echo $jid16 | cut -d' ' -f 4`
cd ../a-1050
jid17=$(sbatch --dependency=afterany:$jid15 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P0.sh)
jid17=`echo $jid17 | cut -d' ' -f 4`
jid18=$(sbatch --dependency=afterany:$jid16 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P1.sh)
jid18=`echo $jid18 | cut -d' ' -f 4`
jid19=$(sbatch --dependency=afterany:$jid17 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P2.sh)
jid19=`echo $jid19 | cut -d' ' -f 4`
jid20=$(sbatch --dependency=afterany:$jid18 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P3.sh)
jid20=`echo $jid20 | cut -d' ' -f 4`
jid21=$(sbatch --dependency=afterany:$jid19 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P4.sh)
jid21=`echo $jid21 | cut -d' ' -f 4`
jid22=$(sbatch --dependency=afterany:$jid20 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P5.sh)
jid22=`echo $jid22 | cut -d' ' -f 4`
jid23=$(sbatch --dependency=afterany:$jid21 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P6.sh)
jid23=`echo $jid23 | cut -d' ' -f 4`
jid24=$(sbatch --dependency=afterany:$jid22 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1050_P7.sh)
jid24=`echo $jid24 | cut -d' ' -f 4`
cd ../a-1055
jid25=$(sbatch --dependency=afterany:$jid23 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P0.sh)
jid25=`echo $jid25 | cut -d' ' -f 4`
jid26=$(sbatch --dependency=afterany:$jid24 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P1.sh)
jid26=`echo $jid26 | cut -d' ' -f 4`
jid27=$(sbatch --dependency=afterany:$jid25 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P2.sh)
jid27=`echo $jid27 | cut -d' ' -f 4`
jid28=$(sbatch --dependency=afterany:$jid26 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P3.sh)
jid28=`echo $jid28 | cut -d' ' -f 4`
jid29=$(sbatch --dependency=afterany:$jid27 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P4.sh)
jid29=`echo $jid29 | cut -d' ' -f 4`
jid30=$(sbatch --dependency=afterany:$jid28 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P5.sh)
jid30=`echo $jid30 | cut -d' ' -f 4`
jid31=$(sbatch --dependency=afterany:$jid29 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P6.sh)
jid31=`echo $jid31 | cut -d' ' -f 4`
jid32=$(sbatch --dependency=afterany:$jid30 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1055_P7.sh)
jid32=`echo $jid32 | cut -d' ' -f 4`
cd ../a-1060
jid33=$(sbatch --dependency=afterany:$jid31 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P0.sh)
jid33=`echo $jid33 | cut -d' ' -f 4`
jid34=$(sbatch --dependency=afterany:$jid32 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P1.sh)
jid34=`echo $jid34 | cut -d' ' -f 4`
jid35=$(sbatch --dependency=afterany:$jid33 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P2.sh)
jid35=`echo $jid35 | cut -d' ' -f 4`
jid36=$(sbatch --dependency=afterany:$jid34 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P3.sh)
jid36=`echo $jid36 | cut -d' ' -f 4`
jid37=$(sbatch --dependency=afterany:$jid35 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P4.sh)
jid37=`echo $jid37 | cut -d' ' -f 4`
jid38=$(sbatch --dependency=afterany:$jid36 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P5.sh)
jid38=`echo $jid38 | cut -d' ' -f 4`
jid39=$(sbatch --dependency=afterany:$jid37 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P6.sh)
jid39=`echo $jid39 | cut -d' ' -f 4`
jid40=$(sbatch --dependency=afterany:$jid38 Slurm_Transfer_08-09-20_rA_AR1C_C_a-1060_P7.sh)
jid40=`echo $jid40 | cut -d' ' -f 4`
