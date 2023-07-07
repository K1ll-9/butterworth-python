/*
 * filter.c
 *
 * Copyright (c) 2018 Disi A
 * 
 * Author: Disi A
 * Email: adis@live.cn
 *  https://www.mathworks.com/matlabcentral/profile/authors/3734620-disi-a
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "filter.h"
#include <stddef.h>

#if DOUBLE_PRECISION
#define COS cos
#define SIN sin
#define TAN tan
#define COSH cosh
#define SINH sinh
#define SQRT sqrt
#define LOG log
#else
#define COS cosf
#define SIN sinf
#define TAN tanf
#define COSH coshf
#define SINH sinhf
#define SQRT sqrtf
#define LOG logf
#endif

BWLowPass* create_bw_low_pass_filter(int order, FTR_PRECISION s, FTR_PRECISION f) {
    BWLowPass* filter = (BWLowPass *) malloc(sizeof(BWLowPass));
    filter -> n = order/2;
    filter -> A = (FTR_PRECISION *)malloc(filter -> n*sizeof(FTR_PRECISION));
    filter -> d1 = (FTR_PRECISION *)malloc(filter -> n*sizeof(FTR_PRECISION));
    filter -> d2 = (FTR_PRECISION *)malloc(filter -> n*sizeof(FTR_PRECISION));
    filter -> w0 = (FTR_PRECISION *)calloc(filter -> n, sizeof(FTR_PRECISION));
    filter -> w1 = (FTR_PRECISION *)calloc(filter -> n, sizeof(FTR_PRECISION));
    filter -> w2 = (FTR_PRECISION *)calloc(filter -> n, sizeof(FTR_PRECISION));
    
    if (filter->d2 == NULL)
    {
        free_bw_low_pass(filter);
        return NULL;
    }

    FTR_PRECISION a = TAN((FTR_PRECISION)(M_PI * f / s));
    FTR_PRECISION a2 = a * a;
    FTR_PRECISION r;
    
    for(int i=0; i < filter -> n; ++i){
        r = SIN((FTR_PRECISION)(M_PI * (2.0 * i + 1.0) / (4.0 * filter->n)));
        s = (FTR_PRECISION) (a2 + 2.0 * a * r + 1.0);
        filter->A[i] = a2 / s;
        filter->d1[i] = (FTR_PRECISION) (2.0 * (1 - a2) / s);
        filter->d2[i] = (FTR_PRECISION)(-(a2 - 2.0 * a * r + 1.0) / s);
    }
    return filter;
}

CHELowPass* create_che_low_pass_filter(int n, FTR_PRECISION epsilon, FTR_PRECISION s, FTR_PRECISION f){
    CHELowPass* filter = (CHELowPass *) malloc(sizeof(CHELowPass));
    filter -> m = n/2;
    filter -> A = (FTR_PRECISION *)malloc(filter -> m*sizeof(FTR_PRECISION));
    filter -> d1 = (FTR_PRECISION *)malloc(filter -> m*sizeof(FTR_PRECISION));
    filter -> d2 = (FTR_PRECISION *)malloc(filter -> m*sizeof(FTR_PRECISION));
    filter -> w0 = (FTR_PRECISION *)calloc(filter -> m, sizeof(FTR_PRECISION));
    filter -> w1 = (FTR_PRECISION *)calloc(filter -> m, sizeof(FTR_PRECISION));
    filter -> w2 = (FTR_PRECISION *)calloc(filter -> m, sizeof(FTR_PRECISION));

    if (filter->d2 == NULL)
    {
        free_che_low_pass(filter);
        return NULL;
    }

    FTR_PRECISION a = TAN((FTR_PRECISION) (M_PI * f/ s));
    FTR_PRECISION a2 = a * a;
    
    FTR_PRECISION u = LOG((FTR_PRECISION) (1.0 + SQRT((FTR_PRECISION) (1.0 + epsilon * epsilon)) / epsilon));
    FTR_PRECISION su = SINH(u/(FTR_PRECISION)n);
    FTR_PRECISION cu = COSH(u/(FTR_PRECISION)n);
    FTR_PRECISION b,c;
    
    int i;
    for(i=0; i<filter->m; ++i){
        b = SIN((FTR_PRECISION) (M_PI * (2.0 * i + 1.0) / (2.0 * n))) * su;
        c = COS((FTR_PRECISION) (M_PI * (2.0 * i + 1.0) / (2.0 * n))) * cu;
        c = b*b + c*c;
        s = (FTR_PRECISION) (a2 * c + 2.0 * a * b + 1.0);
        filter->A[i] = (FTR_PRECISION) (a2 / (4.0 * s));
        filter->d1[i] = (FTR_PRECISION) (2.0 * (1 - a2 * c) / s);
        filter->d2[i] = (FTR_PRECISION) (- (a2 * c - 2.0 * a * b + 1.0) / s);
    }
    filter->ep = (FTR_PRECISION) (2.0 / epsilon);  // used to normalize
    
    return filter;
}


FTR_PRECISION bw_low_pass(BWLowPass* filter, FTR_PRECISION x){
    int i;
    for(i=0; i<filter->n; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i] * (filter->w0[i] + 2.0f * filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x;
}
FTR_PRECISION che_low_pass(CHELowPass* filter, FTR_PRECISION x){
    int i;
    for(i=0; i<filter->m; ++i){
        filter->w0[i] = filter->d1[i]*filter->w1[i] + filter->d2[i]*filter->w2[i] + x;
        x = filter->A[i] * (filter->w0[i] + 2.0f * filter->w1[i] + filter->w2[i]);
        filter->w2[i] = filter->w1[i];
        filter->w1[i] = filter->w0[i];
    }
    return x * filter->ep;
}

FTR_PRECISION softmax(FTR_PRECISION* data, int size, int target_ind){
    FTR_PRECISION sum = 0;
    for(int i = 0; i < size; i++) sum += data[i];
    return data[target_ind]/sum;
}

void free_bw_low_pass(BWLowPass* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter);
}
void free_che_low_pass(CHELowPass* filter){
    free(filter -> A);
    free(filter -> d1);
    free(filter -> d2);
    free(filter -> w0);
    free(filter -> w1);
    free(filter -> w2);
    free(filter);
}