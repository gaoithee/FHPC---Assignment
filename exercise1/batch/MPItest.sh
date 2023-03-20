#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=ex1
#SBATCH --nodes=2
#SBATCH --nodelist=epyc[001-002]
#SBATCH --exclusive
#SBATCH --time=02:00:00

module load architecture/AMD
module load openMPI/4.1.4/gnu/12.2.1

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=1

mpirun -np 60 ./main.x -i -k20000

echo cores,close,k20k,n50,2nodes,no scatter >> mpi20k.dat
for i in {1..256..2}
        do for j in {1..4}
        do mpirun -np $i --map-by core ./main.x -r -n50 -s0 -e1 >> mpi20k.dat
        done
done
