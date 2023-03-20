#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

void write_pgm_image( void *image, int maxval, long xsize, long ysize, const char *image_name);

#define MAXVAL 255

int generate_seed(int omp_rank, int mpi_rank){
  return omp_rank+mpi_rank+omp_rank*mpi_rank;
}


void initialize_serial(const char* filename, unsigned char * world, long size){
    
    //size = number of elements in a row or in a column (a square matrix is used)
    world = (unsigned char *)malloc(size*size*sizeof(unsigned char));

    for(long long i=0; i<size*size; i++){
        
        int val = rand()%100;
        if(val>30){
            world[i]=0; //black=dead
        }else{
            world[i]=MAXVAL; //white=alive
        }
    }

    write_pgm_image(world, MAXVAL, size, size, filename);

    free(world);
}


void initialize_parallel(unsigned char * world, int pRank, unsigned int* rcounts, unsigned int* displs){

   unsigned char * process_world = (unsigned char *)malloc(rcounts[pRank]*sizeof(unsigned char));

  #pragma omp parallel
  {
    //set different seed for each process
    srand(generate_seed(pRank, omp_get_thread_num()));

    #pragma omp for schedule(static)
    for(long long i=0; i<rcounts[pRank]; i++){
        
        int val = rand()%100;
        if(val>30){
            process_world[i]=0; //black = dead
        }else{
            process_world[i]=MAXVAL; //white = alive
        }
    }
    }


MPI_Gatherv(process_world, rcounts[pRank], MPI_UNSIGNED_CHAR, world, rcounts, displs, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

   free(process_world);
}


void choose_initialization(const char * filename, long size, int * argc, char ** argv[]){

    int pRank, pSize; 
unsigned char* world = (unsigned char*)malloc(size*size*sizeof(unsigned char));

    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
    MPI_Comm_size(MPI_COMM_WORLD, &pSize);
    


//here the determination of rcounts and displs:
    unsigned int* displs = (unsigned int *)malloc(pSize*sizeof(unsigned int)); 
    unsigned int* rcounts = (unsigned int *)malloc(pSize*sizeof(unsigned int)); 

long smaller_size;
long cumulative=0;

if(pRank==0){

	for(int i=0; i<pSize; i++){
		//number of rows each process will have to manage
		smaller_size = size%pSize <= i? size/pSize: size/pSize+1;

		//number of elements the process will have to manage
		rcounts[i] = smaller_size*size;

		//index of the first element the process will have to manage
		displs[i] = cumulative;
		
		// tells from which element the next process will have to start
		cumulative = cumulative+rcounts[i];

	}

}

MPI_Bcast(rcounts, pSize, MPI_INT,0, MPI_COMM_WORLD);
MPI_Bcast(displs, pSize, MPI_INT, 0, MPI_COMM_WORLD);

if(pSize > 1){

       	initialize_parallel(world, pRank, rcounts, displs);

	MPI_Barrier(MPI_COMM_WORLD);
	if(pRank==0){
		write_pgm_image(world, MAXVAL, size, size, filename);
	}
}else{

        initialize_serial(filename, world, size);
}
  MPI_Finalize();
  free(world);
}
