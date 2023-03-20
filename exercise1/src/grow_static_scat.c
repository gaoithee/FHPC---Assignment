#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <omp.h>

#if defined(_OPENMP)
#define CPU_TIME ({struct  timespec ts; clock_gettime( CLOCK_REALTIME, &ts ),\
                                          (double)ts.tv_sec +           \
                                          (double)ts.tv_nsec * 1e-9;})
#else

#define CPU_TIME ({struct  timespec ts; clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ),\
                                          (double)ts.tv_sec +           \
                                          (double)ts.tv_nsec * 1e-9;})

#endif

#ifndef RW_PGM
#define RW_PGM
void read_pgm_image( unsigned char  **image, int *maxval, long *xsize, long *ysize, const char *image_name);
void write_pgm_image( unsigned char *image, int maxval, long xsize, long ysize, const char *image_name);
#endif



#define MAXVAL 255


/*
######################################################
#                                                    #
#                                                    #
#               SERIAL FUNCTIONS                     #
#                                                    #
#                                                    #
######################################################
*/

void evaluate_world_serial(unsigned char* world, unsigned char* new_world, long sizex, long sizey){ 
    //funzione che usa world solo tramite lettura, per stabilire se aggiornare o meno new_world
    
    double invsizex=1.0/sizex;
    double invMV=1.0/MAXVAL;
    long square=sizex*sizex;
    for(long k=0; k<sizex*sizey; k++){
        //determine the row and column where we are
        long row = k*invsizex;
        long col = k%sizey;
        //count the number of alive cells in the neighbours
        int SUM = world[(sizex+row-1)%square + (sizey+col-1)%sizey] 
                + world[(sizex+row+0)%square + (sizey+col-1)%sizey] 
                + world[(sizex+row+1)%square + (sizey+col-1)%sizey] 
                + world[(sizex+row-1)%square + (sizey+col+0)%sizey] 
                + world[(sizex+row+1)%square + (sizey+col+0)%sizey] 
                + world[(sizex+row-1)%square + (sizey+col+1)%sizey] 
                + world[(sizex+row+0)%square + (sizey+col+1)%sizey] 
                + world[(sizex+row+1)%square + (sizey+col+1)%sizey];
        SUM = SUM*invMV;
     	//update the world
        new_world[k] = 0;
        if(SUM==2 || SUM==3){
            new_world[k] = MAXVAL;
        }
    }
}


void grw_serial_static(unsigned char* world, long size, int snap, int times){

    //auxiliary world to do the switch
    unsigned char* new_world = (unsigned char *)malloc(size*size*sizeof(unsigned char));

    for(int i=0; i<times; i++){
        unsigned char * ptr1=(i%2==0)? world: new_world; 
        unsigned char * ptr2=(i%2==0)? new_world: world; 

        evaluate_world_alt(ptr1, ptr2, size, size); 
        if(i%snap==0){
            char * fname = (char*)malloc(60);
            sprintf(fname, "snap/snapshot_STATIC_%03d",i+1);
            write_pgm_image(ptr2, MAXVAL, size, size, fname);
            free(fname);
        }    
    }

    free(new_world);
}
//##############################################################################################################################

/*
######################################################
#                                                    #
#                                                    #
#                OTHER FUNCTIONS                     #
#                                                    #
#                                                    #
######################################################
unsigned char checksingle(unsigned char *matrix,  long i, long N,long M){
    int left=N*(i%N==0);
    int right=N*((i+1)%N==0);
    long square=N*M;
    long up=square*((i-M)<0);
    long down=square*((i+M)>=square);
    unsigned char status=(
        matrix[i-1 + left]  // west el
        +matrix[i+1 - right] // east el
        +matrix[i-M + up] //north el
        +matrix[i+M - down] // south el
        +matrix[i-1-M + left + up] // nort-west el
        +matrix[i+1-M - right + up] // north-east el
        +matrix[i-1+M + left - down] // south-west el
        +matrix[i+1+M - right - down] // south-east el
        )/255; // note: all the elements need a correction to account for the fact that we want a 
              // "pacman effect": left and right apply a correction on the elements of the first and last 
              // column, up and down apply a correction on the elements of the first and last row
    return ((status==5 || status==6)? 0: MAXVAL );
}
void updatematr(unsigned char *world,unsigned char *new_world,long sizex,long sizey){ //generalized for single pieces
        
        for(long i=0;i<sizex*sizey;i++){
            new_world[i]=checksingle(world, i, sizex,sizey);
        }
    
} 
void printmatrix(unsigned char* matrix, int N,int M){
    for (int i=0;i<N;i++){
        for (int j=0;j<M;j++){
            printf("%d ",matrix[i*N+j]);
        }
        printf("\n");
    }
    printf("\n");
} 
*/

void evaluate_world_alt(unsigned char* world, unsigned char* new_world, long sizex, long sizey){ 

    double invsizex=1.0/sizex;
    double invMV=1.0/MAXVAL;
	#pragma omp parallel for
    for(long k=sizex; k<sizex*(sizey-1); k++){

        long col = k%sizex;
        long r = k*invsizex;
    
        //Calculate the neighbours
        long col_prev = col-1>=0 ? col-1 : sizex-1;
        long col_next = col+1<sizex ? col+1 : 0;
        long r_prev = r-1;
        long r_next = r+1;

        //Determine the number of alive neighbours
        int sum = world[r_prev*sizex+col_prev]+
                  world[r_prev*sizex+col]+
                  world[r_prev*sizex+col_next]+
                  world[r*sizex+col_prev]+
                  world[r*sizex+col_next]+
                  world[r_next*sizex+col_prev]+
                  world[r_next*sizex+col]+
                  world[r_next*sizex+col_next];
        sum = sum*invMV;

        //Update the cell
        new_world[k] = 0;
        if((sum==2) || (sum==3)){
            new_world[k]=MAXVAL;
        }
    
    }
} 

//###########################################################################################################################



/*
######################################################
#                                                    #
#                                                    #
#               PARALLEL FUNCTIONS                   #
#                                                    #
#                                                    #
######################################################
*/

void evaluate_world(unsigned char* world, unsigned char* new_world, long sizex, long sizey, int times, int pRank, int pSize, MPI_Status* status, MPI_Request* req){ 

    double invsizex=1.0/sizex;
    double invMV=1.0/MAXVAL;
    #pragma omp parallel for
    for(long long k=sizex; k<sizex*(sizey-1); k++){

        long col = k%sizex;
        long r = k*invsizex;
    
        //Calculate the neighbours
        long col_prev = col-1>=0 ? col-1 : sizex-1;
        long col_next = col+1<sizex ? col+1 : 0;
        long r_prev = r-1;
        long r_next = r+1;

        //Determine the number of dead neighbours
        int sum = world[r_prev*sizex+col_prev]+
                  world[r_prev*sizex+col]+
                  world[r_prev*sizex+col_next]+
                  world[r*sizex+col_prev]+
                  world[r*sizex+col_next]+
                  world[r_next*sizex+col_prev]+
                  world[r_next*sizex+col]+
                  world[r_next*sizex+col_next];
        sum = sum*invMV;
    
        //Update the cell
        new_world[k] = 0;
        if((sum==2) || (sum==3)){
            new_world[k]=MAXVAL;
        }
    
    }
    
    // now let's send the needed rows
    int tag_odd = 2*times;
    int tag_even = 2*times+1;
    MPI_Barrier(MPI_COMM_WORLD);
    //each process sends his first and last row to respectively the process with rank-1 and rank + 1.
    //Process 0 send his fist line to process size -1 
    //Process size-1 send his last row to process 0
    #pragma omp master
    {
    if(pRank == pSize-1){

        MPI_Isend(&new_world[(sizey-2)*sizex], sizex, MPI_UNSIGNED_CHAR, 0, tag_even, MPI_COMM_WORLD, req);
        MPI_Isend(&new_world[sizex], sizex, MPI_UNSIGNED_CHAR, pRank-1, tag_odd, MPI_COMM_WORLD, req);
      
        MPI_Recv(new_world, sizex, MPI_UNSIGNED_CHAR, pRank-1, tag_even, MPI_COMM_WORLD, status);
        MPI_Recv(&new_world[(sizey-1)*sizex], sizex, MPI_UNSIGNED_CHAR, 0, tag_odd, MPI_COMM_WORLD, status);
    }


    if(pRank == 0){
    
        MPI_Isend(&new_world[(sizey-2)*sizex], sizex, MPI_UNSIGNED_CHAR, 1, tag_even, MPI_COMM_WORLD, req);
        MPI_Isend(&new_world[sizex], sizex, MPI_UNSIGNED_CHAR, pSize-1, tag_odd, MPI_COMM_WORLD, req);

        MPI_Recv(new_world, sizex, MPI_UNSIGNED_CHAR, pSize-1, tag_even, MPI_COMM_WORLD, status);
        MPI_Recv(&new_world[(sizey-1)*sizex], sizex, MPI_UNSIGNED_CHAR, 1, tag_odd, MPI_COMM_WORLD, status);
    }

    if((pRank != 0) & (pRank != pSize-1)){

        MPI_Isend(&new_world[(sizey-2)*sizex], sizex, MPI_UNSIGNED_CHAR, pRank+1, tag_even, MPI_COMM_WORLD, req);
        MPI_Isend(&new_world[sizex], sizex, MPI_UNSIGNED_CHAR, pRank-1, tag_odd, MPI_COMM_WORLD, req);

        MPI_Recv(&new_world[(sizey-1)*sizex], sizex, MPI_UNSIGNED_CHAR, pRank+1, tag_odd, MPI_COMM_WORLD, status);
        MPI_Recv(new_world, sizex, MPI_UNSIGNED_CHAR, pRank-1, tag_even, MPI_COMM_WORLD, status);
    }
    }
	
}


void grw_parallel_static(unsigned char* world, long size, int pSize, int pRank, unsigned int* scounts, unsigned int* displs, unsigned int* rcounts_g, unsigned int* displs_g, int snap, int times){
    MPI_Status status;
    MPI_Request req;
    
    //auxiliary local worlds for each process
    //new_world collects the output of Scatterv, temp_new_world is the second world needed to do the switch
    unsigned char* new_world = (unsigned char *)malloc(scounts[pRank]*sizeof(unsigned char));
	unsigned char* temp_new_world = (unsigned char *)malloc(scounts[pRank]*sizeof(unsigned char));
	double tstart = CPU_TIME;    
    MPI_Scatterv(world, scounts, displs, MPI_UNSIGNED_CHAR, new_world, scounts[pRank], MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    	double tend = CPU_TIME;
    if(pRank==0){
	printf("%f ",tend-tstart);  //check time needed by Scatterv call
    }
    double invsize=1.0/size;
    for(int i=0; i<times; i++){

        unsigned char * ptr1=(i%2==0)? new_world: temp_new_world;
        unsigned char * ptr2=(i%2==0)? temp_new_world: new_world;

        evaluate_world(ptr1, ptr2, size, (long)scounts[pRank]*invsize, times, pRank, pSize, &status, &req); 

        if(i%snap==0){
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gatherv(&ptr2[size], rcounts_g[pRank], MPI_UNSIGNED_CHAR, world, rcounts_g, displs_g, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

            if(pRank==0){
                char * fname = (char*)malloc(60);
                sprintf(fname, "snap/snapshot_STATIC_%03d",i+1);
	            write_pgm_image(world, MAXVAL, size, size, fname);
                free(fname);
            }

        }
    }

    free(new_world);
    free(temp_new_world);
}
//###############################################################################################################################

/*
######################################################
#                                                    #
#                                                    #
#           MAIN FUNCTION OF THE PROGRAM             #
#                                                    #
#                                                    #
######################################################
*/
void run_static(char * filename, int times, int dump, int * argc, char ** argv[]){

    unsigned char* world;
    long size=0;
    int maxval=0;

    int pRank, pSize;
    int mpi_provided_thread_level;
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED,&mpi_provided_thread_level);
    if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) {
	    printf("a problem arose when asking for MPI_THREAD_FUNNELED level\n");
	    MPI_Finalize();
	    exit( 1 );
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
    MPI_Comm_size(MPI_COMM_WORLD, &pSize);  

	double tstart = CPU_TIME;

    //serial read
    read_pgm_image(&world, &maxval, &size, &size, filename);
    if (pSize>size) {
        if(pRank==0){
            printf("Program interrupted. You ran with %d processes on a matrix of size %ld x %ld, set the number of processes to a value lower or equal than %ld \n",pSize,size,size,size);
        }
        exit(1);
    }
    //temporary world that contains two more lines, needed for the parallel case
    unsigned char* temp_world = (unsigned char *)malloc(size*(size+2)*sizeof(unsigned char));
	if(pRank==0){

    //let's fill the world
    #pragma omp parallel for
    for(long i=0; i<size*(size+2); i++){
    
        if((i>= size) & (i<size*(size+1))){
            temp_world[i] = world[i-size];
    
        }else if(i<size){
	        temp_world[i]=world[size*(size-1)+i];
        
        }else{
	        temp_world[i]=world[i-size*(size+1)];
        }
    
    }
}

    //auxiliary vectors for Scatterv
    unsigned int* displs = (unsigned int *)malloc(pSize*sizeof(unsigned int));  //starting index for each process
    unsigned int* scounts = (unsigned int *)malloc(pSize*sizeof(unsigned int)); //number of elements to assign to each process

    //auxiliary vectors for Gatherv
    unsigned int* displs_g = (unsigned int *)malloc(pSize*sizeof(unsigned int));  //starting index for each process
    unsigned int* rcounts_g = (unsigned int *)malloc(pSize*sizeof(unsigned int)); //number of elements to assign to each process


    if(pRank==0){  //no need to repeat this for all the processes
        unsigned int smaller_size;
        unsigned int cumulative=0;
        unsigned int cumulative_g=0;
        unsigned int std_size=size/pSize;

	    for(int i=0; i<pSize; i++){
		
            smaller_size = size%pSize <= i? std_size: std_size+1; //work for each process
            // also smaller_size= size/pSize+(size%pSize>pRank)

		    scounts[i]=(smaller_size+2)*size;
		    displs[i]=cumulative;
            cumulative = cumulative+(scounts[i]-2*size);
            
            rcounts_g[i]=smaller_size*size;
		    displs_g[i]=cumulative_g;
		    cumulative_g = cumulative_g+rcounts_g[i];
		   

    	}

    }

    // now let's spread the vectors we just created
    MPI_Bcast(scounts, pSize, MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(displs, pSize, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(rcounts_g, pSize, MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(displs_g, pSize, MPI_INT, 0, MPI_COMM_WORLD);

    
    // main call
    if(pSize > 1){
	    grw_parallel_static(temp_world, size, pSize, pRank, scounts, displs, rcounts_g, displs_g, dump, times); 
    
    }else{
        grw_serial_static(world, size, dump, times);
    }
    

    free(temp_world);
    free(world);
    double tend = CPU_TIME;
	if(pRank==0){
    printf("%d %f\n",pSize,tend-tstart);  //total time needed
}    
MPI_Finalize();
}
