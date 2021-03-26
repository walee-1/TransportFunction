#!/bin/bash
jid1=$(sbatch Slurm_Transfer_08-09-20_NC_opt_d0_C_b15_P0.sh)
jid2=$(sbatch Slurm_Transfer_08-09-20_NC_opt_d0_C_b15_P1.sh)
jid1=`echo $jid1 | cut -d' ' -f 4`
jid2=`echo $jid2 | cut -d' ' -f 4`
jid3=$(sbatch --dependency=afterany:$jid1 Slurm_Transfer_08-09-20_NC_opt_d0_C_b15_P2.sh)
jid3=`echo $jid3 | cut -d' ' -f 4`
jid4=$(sbatch --dependency=afterany:$jid2 Slurm_Transfer_08-09-20_NC_opt_d0_C_b15_P3.sh)
jid5=$(sbatch --dependency=afterany:$jid3 Slurm_Transfer_08-09-20_NC_opt_d0_C_b15_P4.sh)
