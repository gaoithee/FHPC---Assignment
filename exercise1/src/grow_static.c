#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
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

void evaluate_world_alt(unsigned char* world, unsigned char* new_world, long sizex, long sizey){ 

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
} 

void grw_serial_static(unsigned char* world, long size, int snap, int times){

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
    return ((status==2 || status==3)? MAXVAL : 0);
}
void updatematr(unsigned char *world,unsigned char *new_world,long sizex,long sizey){ //generalized for single pieces
        
        for(long i=0;i<sizex*sizey;i++){
            new_world[i]=checksingle(world, i, sizex,sizey);
        }    
} 

void printmatrix(unsigned char* matrix, int N,int M, int prank){
    printf("p%d\n",prank);
    for (int i=0;i<N;i++){
        for (int j=0;j<M;j++){
            printf("%d ",matrix[i*N+j]);
        }
        printf("\n");
    }
    printf("\n");
} 


void evaluate_world_serial(unsigned char* world, unsigned char* new_world, long sizex, long sizey){ 
    //funzione che usa world solo tramite lettura, per stabilire se aggiornare o meno new_world
    
    double invsizex=1.0/sizex;
    double invMV=1.0/MAXVAL;
    for(long k=0; k<sizex*sizey; k++){
        //determiniamo il numero di riga e di colonna in cui siamo
        long row = k*invsizex;
        long col = k%sizey;
        long square=sizex*sizex;
        //contiamo il numero di vive e di morte nelle circostanze
        int SUM = world[(sizex+row-1)%square + (sizey+col-1)%sizey] 
                + world[(sizex+row+0)%square + (sizey+col-1)%sizey] 
                + world[(sizex+row+1)%square + (sizey+col-1)%sizey] 
                + world[(sizex+row-1)%square + (sizey+col+0)%sizey] 
                + world[(sizex+row+1)%square + (sizey+col+0)%sizey] 
                + world[(sizex+row-1)%square + (sizey+col+1)%sizey] 
                + world[(sizex+row+0)%square + (sizey+col+1)%sizey] 
                + world[(sizex+row+1)%square + (sizey+col+1)%sizey];
        //dividendo per il valore (MAXVAL=255) che rappresenta le celle morte
        //otteniamo il numero di morti che circondano la cella in posizione (row, col)
        SUM = SUM*invMV;
        //di default, la facciamo morire
        new_world[k] = 0;
        //ma se per caso è circondata da 2 o 3 celle vive sulle 8 vicine,
        //che è analogo a dire che è circondata da 5 o 6 celle morte, allora vive
        if(SUM==2 || SUM==3){
            new_world[k] = MAXVAL;
        }
    }
}
*/

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


void grw_parallel_static(unsigned char* process_world, unsigned char* world, long size, int pSize, int pRank, unsigned int* rcounts_g, unsigned int* displs_g, int snap, int times){
    MPI_Status status;
    MPI_Request req;
    
    unsigned char* new_world = (unsigned char *)malloc((rcounts_g[pRank]+2*size)*sizeof(unsigned char));
    double invsize=1.0/size;
    for(int i=0; i<times; i++){

        unsigned char * ptr1=(i%2==0)? process_world: new_world; 
        unsigned char * ptr2=(i%2==0)? new_world: process_world; 

        evaluate_world(ptr1, ptr2, size, (long)rcounts_g[pRank]*invsize+2, times, pRank, pSize, &status, &req); 

        if(i%snap==0){
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

    //auxiliary vectors for Gatherv
    unsigned int* displs_g = (unsigned int *)malloc(pSize*sizeof(unsigned int));  //starting index for each process
    unsigned int* rcounts_g = (unsigned int *)malloc(pSize*sizeof(unsigned int)); //number of elements to assign to each process


    if(pRank==0){  //no need to repeat this for all the processes
        long smaller_size;
        long cumulative_g=0;
        long std_size=size/pSize;

	    for(int i=0; i<pSize; i++){
		
            smaller_size = size%pSize <= i? std_size: std_size+1; //work for each process
            // also smaller_size= std_size+(size%pSize>pRank)
            rcounts_g[i]=smaller_size*size;
		    displs_g[i]=cumulative_g;
		    cumulative_g = cumulative_g+rcounts_g[i];
    	}
    }

    MPI_Bcast(rcounts_g, pSize, MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(displs_g, pSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    long smaller_size = rcounts_g[pRank]/size;
    long long beginning = displs_g[pRank] - size;
    
    unsigned char* process_world = (unsigned char*)malloc(size*(smaller_size+2)*sizeof(unsigned char));
    if((pRank!=0) && (pRank!=pSize-1)){
        #pragma omp parallel for
	for(long long k=0; k<size*(smaller_size+2); k++){
            process_world[k] = world[beginning+k];
        }
    }
    if(pRank == 0){
	#pragma omp parallel for
        for(long long k=0; k<size*(smaller_size+2); k++){
                process_world[k] = k<size ? world[(size-1)*size+k] : world[k-size];
         }
    }
    if(pRank == pSize-1){
	#pragma omp parallel for
        for(long long k=0; k<size*(smaller_size+2); k++){
                process_world[k] = k<(size*(smaller_size+1)) ? world[beginning+k] : world[k-size*(smaller_size+1)];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double tstart2=CPU_TIME;    
    // main call
    if(pSize > 1){
	    grw_parallel_static(process_world, world, size, pSize, pRank, rcounts_g, displs_g, dump, times); 
    
    }else{
        grw_serial_static(world, size, dump, times);
    }
    

    free(process_world);
    free(world);
    free(displs_g);
    free(rcounts_g);
    double tend = CPU_TIME;
	if(pRank==0){
    printf("%d %f %f\n",pSize,tend-tstart2,tend-tstart);
}    
MPI_Finalize();
}

