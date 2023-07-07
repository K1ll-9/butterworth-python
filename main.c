#include <stdio.h>
#include <stdlib.h>
#include "filter.h"

#define DATA_LEN  2205
#define RUNS      5
#define ORDER     4


int main() {

  BWLowPass* bpfilter;

  double      SR=44100;
  double      filtered_signal[RUNS][DATA_LEN];
  double      input[DATA_LEN];

  printf("Opening and reading data file...\n");

  FILE *fp = fopen("inputs.raw", "r");
  if(fp == NULL) {
    perror("Unable to open file!");
    exit(1);
  }
  char line[30];

  //double filter_lowfc=60;
  double filter_highfc=205;

  bpfilter = create_bw_low_pass_filter(ORDER, SR, filter_highfc);

  unsigned int i=0;
  while(fgets(line, sizeof(line), fp) != NULL) {
    sscanf(line, "%lf", &input[i]);
    i++;
  }

  fclose(fp);

  printf("Data stored into input[].\n");

  for(int j=0; j<RUNS; j++) {
    printf("Run #%d:\n", RUNS);
    for(int i=0; i<DATA_LEN; i++) {
      filtered_signal[j][i] = bw_low_pass(bpfilter, input[i]);
    }
  }


  printf("Comparing results:\n");
  for(int i=0; i<500; i++) {
    printf("[bp_filter %i] input: %lf | ", i, input[i]);
    for(int j=0; j<RUNS; j++)
      printf(" f_out#%d %lf |", j, filtered_signal[j][i]);
    printf("\n");
  }

/*
    filter -> n = order/2;
    filter -> A = (FTR_PRECISION *)malloc(filter -> n*sizeof(FTR_PRECISION));
    filter -> d1 = (FTR_PRECISION *)malloc(filter -> n*sizeof(FTR_PRECISION));
    filter -> d2 = (FTR_PRECISION *)malloc(filter -> n*sizeof(FTR_PRECISION));
    filter -> w0 = (FTR_PRECISION *)calloc(filter -> n, sizeof(FTR_PRECISION));
    filter -> w1 = (FTR_PRECISION *)calloc(filter -> n, sizeof(FTR_PRECISION));
    filter -> w2 = (FTR_PRECISION *)calloc(filter -> n, sizeof(FTR_PRECISION));
*/
  printf("\nFilter settings:\n");
  printf("Order: %d\n", ORDER);
  printf("n: %d\n", bpfilter->n);
  printf("d1: %lf\n", *bpfilter->d1);
  printf("d2: %lf\n", *bpfilter->d2);
  printf("w0[0]: %lf %lf\n", bpfilter->w0[0], bpfilter->w0[1]);
  printf("w1[0]: %lf %lf\n", bpfilter->w1[0], bpfilter->w1[1]);
  printf("w2[0]: %lf %lf\n", bpfilter->w2[0], bpfilter->w2[1]);


  free_bw_low_pass(bpfilter);

  return(0);
}
