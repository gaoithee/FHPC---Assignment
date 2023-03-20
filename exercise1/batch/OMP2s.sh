#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=ex1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00

module load architecture/AMD
module load openMPI/4.1.4/gnu/12.2.1

export OMP_PLACES=cores
export OMP_PROC_BIND=close

mpirun -np 2 ./main.x -i -k25000

echo cores,close,k20k,n50,1node >> omp2s25k.dat
for i in {1..64}
	do
	 export OMP_NUM_THREADS=$i
	for j in {1..5}
	do mpirun -np 2 --map-by socket ./main.x -r -n50 -s0 -e1 >> omp2s25k.dat
	done
done



