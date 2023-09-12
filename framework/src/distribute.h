#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H


int distribute_init(int num_columns, int *myrank, int *col_s, int *col_e);

int distribute_final(void);


#endif
