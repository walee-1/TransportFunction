#!/bin/bash
jid1=$(sbatch Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P0.sh)
jid1=`echo $jid1 | cut -d' ' -f 4`
jid2=$(sbatch Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P1.sh)
jid2=`echo $jid2 | cut -d' ' -f 4`
jid3=$(sbatch --dependency=afterany:$jid1 Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P2.sh)
jid3=`echo $jid3 | cut -d' ' -f 4`
jid4=$(sbatch --dependency=afterany:$jid2 Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P3.sh)
jid4=`echo $jid4 | cut -d' ' -f 4`
jid5=$(sbatch --dependency=afterany:$jid3 Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P4.sh)
jid5=`echo $jid5 | cut -d' ' -f 4`
jid6=$(sbatch --dependency=afterany:$jid4 Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P5.sh)
jid6=`echo $jid6 | cut -d' ' -f 4`
jid7=$(sbatch --dependency=afterany:$jid5 Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P6.sh)
jid7=`echo $jid7 | cut -d' ' -f 4`
jid8=$(sbatch --dependency=afterany:$jid6 Slurm_Transfer_08-09-20_SC_Prot_opt_drF_C_a-1060_P7.sh)
