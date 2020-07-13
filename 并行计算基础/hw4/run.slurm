#!/bin/bash

#SBATCH -o Hw.out
#SBATCH --partition=gpu
#SBATCH -J mat_mul 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 10:00
m = 100, n = 101， s = 311;
./mat_mul $m $n $s
