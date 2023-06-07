#!/bin/bash
#
#
#SBATCH --job-name=Mma_TransP1
#SBATCH --output=MmaOut_%A.txt
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=medium
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=3500mb
#

BASEPATH="/users/waleed.khalid/Mma/MergerTransport/"
PLACEHOLDERPATH=NC_Prot_opt_drA/a-1050/

date
module load mathematica/12.1.0
cd $BASEPATH$PLACEHOLDERPATH
math -script Transfer_08-09-20_NC_Prot_opt_drA_C_a-1050_P3.m
date