#ifndef included_mpi_cplx_hpp
#define included_mpi_cplx_hpp

#include "mpi.h"
#include "complexify.h"
  
int MPI_Send_Cplx( void *buf_re, void *buf_im, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm );  

int MPI_Irecv_Cplx( void *buf_re, void *buf_im, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm,  MPI_Request *request );  

#endif

