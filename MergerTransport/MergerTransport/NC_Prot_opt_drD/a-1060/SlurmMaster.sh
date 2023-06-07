#!/bin/bash
jid1=$(sbatch Slurm_Transfer_08-09-20_NC_Prot_opt_drD_C_a-1060_P3.sh)
jid1=`echo $jid1 | cut -d' ' -f 4`
jid2=$(sbatch Slurm_Transfer_08-09-20_NC_Prot_opt_drD_C_a-1060_P4.sh)
jid2=`echo $jid2 | cut -d' ' -f 4`
jid3=$(sbatch --dependency=afterany:$jid1 Slurm_Transfer_08-09-20_NC_Prot_opt_drD_C_a-1060_P5.sh)
jid3=`echo $jid3 | cut -d' ' -f 4`
jid4=$(sbatch --dependency=afterany:$jid2 Slurm_Transfer_08-09-20_NC_Prot_opt_drD_C_a-1060_P6.sh)
jid4=`echo $jid4 | cut -d' ' -f 4`
jid5=$(sbatch --dependency=afterany:$jid3 Slurm_Transfer_08-09-20_NC_Prot_opt_drD_C_a-1060_P7.sh)
