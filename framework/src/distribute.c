#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"


#define mpi_check(err) \
{ \
    if (err != MPI_SUCCESS) \
    { \
        int buffer_len; \
        char buffer[128]; \
        MPI_Error_string(err, buffer, &buffer_len); \
        fprintf(stderr, "MPI error: %s\n", buffer); \
        return 2; \
    } \
}


int distribute_init(int num_columns, int *myrank, int *col_s, int *col_e)
{
    /*Initialize MPI.*/
    mpi_check(MPI_Init(NULL, NULL));

    /*Get rank information.*/
    int rank;
    mpi_check(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    *myrank = rank;
    int num_ranks;
    mpi_check(MPI_Comm_size(MPI_COMM_WORLD, &num_ranks));

    /*Calculating starting and ending column indices.*/
    float x = ((float)num_columns)/((float)num_ranks);
    int ndiv = (int)ceil(x);
    int dx = (int)floor(x);
    *col_s = rank*dx;
    *col_e = *col_s + dx;
    if (*col_e >= num_columns)
    {
        *col_e = num_columns - 1;
    }
    return 0;
}


int distribute_final(void)
{
    /*Finalize MPI.*/
    mpi_check(MPI_Finalize());
    return 0;
}
