#!/bin/bash
#
#
#SBATCH --job-name=Mma_TransP1
#SBATCH --output=MmaOut_%A.txt
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=3500mb
#

BASEPATH="/users/daniel.moser/Mma/TransferProject/TransportFunction/MergerTransport/"
PLACEHOLDERPATH='b15/'

date
module load mathematica/12.1.0
cd $BASEPATH$PLACEHOLDERPATH
math -script Transfer_28-08-20_NC_opt_dAlpha_C_b15_P1.m
date
