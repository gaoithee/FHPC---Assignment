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

mpirun -np 1 ./main.x -i -k15000

echo cores,close,k15k,n10,updatematr >> compare.dat
for j in {1..10}
        do mpirun -np 1 ./main.x -r -n10 -s0 -e1 >> compare.dat
done

echo cores,close,k15k,n10,evaluate_world_serial >> compare.dat
for j in {1..10}
	do mpirun -np 1 ./main2.x -r -n10 -s0 -e1 >> compare.dat
done

echo cores,close,k15k,n10,evaluate_world_alt >> compare.dat
for j in {1..10}
        do mpirun -np 1 ./main3.x -r -n10 -s0 -e1 >> compare.dat
done

