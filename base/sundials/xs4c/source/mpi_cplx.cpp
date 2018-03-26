#include "mpi_cplx.h"

int MPI_Send_Cplx( void *buf_re, void *buf_im, int count, MPI_Datatype datatype, int dest,
                   int tag, MPI_Comm comm )
{
  int retval;

  retval = MPI_Send( buf_re, count, datatype, dest, tag, comm );
  if (retval != MPI_SUCCESS) return(retval);

  retval = MPI_Send( buf_im, count, datatype, dest, tag, comm );

  return (retval);
}

int MPI_Irecv_Cplx( void *buf_re, void *buf_im, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm,  MPI_Request *request )
{
  int retval;

  retval =  MPI_Irecv( buf_re, count, datatype, dest, tag, comm, request );
  if (retval != MPI_SUCCESS) return(retval);

  retval = MPI_Irecv( buf_im, count, datatype, dest, tag, comm, request );
  return(retval);
}
