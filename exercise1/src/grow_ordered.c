#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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

void update_world(unsigned char * world, long size){

    // I want to update the matrix in a ordered way:
    //starting from the first element, neighbours will be considered in order to establish whether the cell must be remain alive or not.
    //then i will proceed in row-major direction, just as an old writing machine
    int SUM;
    double invsize=1.0/size;
    double invMV=1.0/MAXVAL;
    for(long long k=0; k<size*size; k++){

      long row = k*invsize;
      long col = k%size;

      SUM = 0;

        
        SUM = world[(size+row-1)%size*size + (size+col-1)%size] 
            + world[(size+row+0)%size*size + (size+col-1)%size] 
            + world[(size+row+1)%size*size + (size+col-1)%size] 
            + world[(size+row-1)%size*size + (size+col+0)%size] 
            + world[(size+row+1)%size*size + (size+col+0)%size] 
            + world[(size+row-1)%size*size + (size+col+1)%size] 
            + world[(size+row+0)%size*size + (size+col+1)%size] 
            + world[(size+row+1)%size*size + (size+col+1)%size];
            

        SUM = SUM*invMV; //number of alive elements that surround the selected cell


        //dead if 1 or 0 cells alive in its neighbours (due to undercrowding) or if 
        //more than 3 cells alive in its neighbours (due to overcrowding)
        world[k] = 0;

        //alive (encoded as 255) if 2 or 3 cells in its neighbours (made of 8 cells) are alive
        if(SUM==2 || SUM==3){
           world[k] = MAXVAL; 
        }

    }
}

/*
void update_world(unsigned char* world, long size){ 

    double invsize=1.0/size;
    double invMV=1.0/MAXVAL;
    for(long long k=0; k<size*size; k++){

        long col = k%size;
        long r = k*invsize;
    
    	//Calculate the neighbours
        long col_prev = col-1>=0 ? col-1 : size-1;
        long col_next = col+1<size ? col+1 : 0;
        long r_prev = r-1;
        long r_next = r+1;

        //Determine the number of dead neighbours
        int sum = world[r_prev*size+col_prev]+
                  world[r_prev*size+col]+
                  world[r_prev*size+col_next]+
                  world[r*size+col_prev]+
                  world[r*size+col_next]+
                  world[r_next*size+col_prev]+
                  world[r_next*size+col]+
                  world[r_next*size+col_next];
        sum = sum*invMV;
    
        //Update the cell
        world[k] = 0;
        if((sum==2) || (sum==3)){
            world[k]=MAXVAL;
        }
    }
} 

*/

//the program asks as input the number of steps that we need to perform the game of life: here "times"
void grw_serial(unsigned char* world, long size, int times, int snap){

  for(int i=1; i<=times; i++){

    //at each iteration, we update the world
    update_world(world, size);

    //snaps = every how many steps a dump of the system is saved on a file (0 meaning only at the end)
    if(i%snap==0){
      //the number of iterations to save in fname could not exceed times
      char * fname = (char*)malloc(60);
      sprintf(fname, "snap/snapshot_ORDERED_%03d",i);      
      write_pgm_image(world, MAXVAL, size, size, fname);
      free(fname);
    }
  }
}


void run_ordered(char * filename, int times, int dump){
  
  double tstart = CPU_TIME;

  unsigned char * world;
  long size = 0;
  int maxval = 0; 
  
  //first of all, i need to read the previous world state
  read_pgm_image(&world, &maxval, &size, &size, filename);
  //main call
  grw_serial(world, size, times, dump);

  double tend = CPU_TIME;
  printf("%f\n",tend-tstart);
  free(world);
}

