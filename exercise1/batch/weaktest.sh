#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=ex1
#SBATCH --nodes=3
#SBATCH --exclusive
#SBATCH --time=02:00:00

module load architecture/AMD
module load openMPI/4.1.4/gnu/12.2.1

export OMP_PLACES=cores
export OMP_PROC_BIND=close


echo cores,close,kvar,n50,no scatter >> weaknoscatt.dat


export OMP_NUM_THREADS=64
mpirun -np 60 main.x -i -k 10000 -f init2
for j in {1..5}
        do mpirun -np 1 --map-by socket main.x -r -n50 -s0 -e1 -f init2 >> weaknoscatt.dat
done

export OMP_NUM_THREADS=64   
mpirun -np 60 main.x -i -k 14143 -f init2
for j in {1..5}
        do mpirun -np 2 --map-by socket main.x -r -n50 -s0 -e1 -f init2 >> weaknoscatt.dat
done

export OMP_NUM_THREADS=64     
mpirun -np 60 main.x -i -k 17321 -f init2
for j in {1..5}
        do mpirun -np 3 --map-by node --bind-to socket main.x -r -n50 -s0 -e1 -f init2 >> weaknoscatt.dat
done

export OMP_NUM_THREADS=64
mpirun -np 60 main.x -i -k 20000 -f init2
for j in {1..5}
        do mpirun -np 4 --map-by node --bind-to socket main.x -r -n50 -s0 -e1 -f init2 >> weaknoscatt.dat
done

export OMP_NUM_THREADS=64
mpirun -np 60 main.x -i -k 22361 -f init2
for j in {1..5}
        do mpirun -np 5 --map-by node --bind-to socket main.x -r -n50 -s0 -e1 -f init2 >> weaknoscatt.dat
done

export OMP_NUM_THREADS=64    
mpirun -np 60 main.x -i -k 24494 -f init2
for j in {1..5}
        do mpirun -np 6 --map-by node --bind-to socket main.x -r -n50 -s0 -e1 -f init2 >> weaknoscatt.dat
done
