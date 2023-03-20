#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=p2-sp-mkl-spr
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00

module load architecture/AMD
module load mkl

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

for numcores in {1..64}
	do export OMP_NUM_THREADS=$numcores
	echo $numcores >> weak-sp-mkl-spr.dat
	for i in {1..10}
		do ./gemm_mkl_sp.x 12000 12000 12000 | grep GFLOPS >> weak-sp-mkl-spr.dat
	done
done
