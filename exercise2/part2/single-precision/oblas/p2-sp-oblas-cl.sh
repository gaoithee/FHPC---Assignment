#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=p2-sp-oblas-cl
#SBATCH --nodes=1
#SBATCH --nodelist=epyc[008]
#SBATCH --exclusive
#SBATCH --time=02:00:00

module load architecture/AMD
module load openBLAS/0.3.21-omp

export OMP_PLACES=cores
export OMP_PROC_BIND=close

for numcores in {1..64}
        do export OMP_NUM_THREADS=$numcores
        echo $numcores >> weak-sp-oblas-cl.dat
        for i in {1..10}
                do ./gemm_oblas_sp.x 12000 12000 12000 | grep GFLOPS >> weak-sp-oblas-cl.dat
        done
done
