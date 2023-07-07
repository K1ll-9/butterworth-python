#include <stdio.h>
#include <math.h>

#define DEBUG  0

typedef struct filter {
  double a0;
  double a1;
  double a2;
  double b0;
  double b1;
  double b2;

  double x0;
  double x1;
  double x2;
  double y1;
  double y2;
} ChebFilter;

ChebFilter call_205(int P, ChebFilter filter, double FC, int NP, int LH, double PR) { 

  double rp= -cos(M_PI/(NP*2)) + (P-1)*M_PI/NP;
  double ip=  sin(M_PI/(NP*2)) + (P-1)*M_PI/NP;

  if(DEBUG) {
    printf("\n[call_205 #%d] rp= %.10lf\n", P, rp);
    printf("[call_205 #%d] ip= %.10lf\n", P, ip);
    printf("\n");
  }

  if(PR!=0) { 
   double es= sqrt((double)100/pow(100-PR,2)-1);
   double vx= (double)(1/NP) * log( (double)(1/es) / sqrt( (double)(1/pow(es,2))+1));
   double kx= (double)(1/NP) * log( (double)(1/es) / sqrt( (double)(1/pow(es,2))-1));
          kx= (exp(kx) + exp(-kx))/2;
          rp= rp * ( (exp(vx) - exp(-vx) ) /2 ) / kx;
          ip= ip * ( (exp(vx) + exp(-vx) ) /2 ) / kx;
  }

  double t=2*tan((double)1/2);
  double w= 2*M_PI*FC;
  double m=pow(rp,2) + pow(ip,2);
  double d= 4 - 4*rp*t + m*pow(t,2);
  double x0= pow(t,2)/d;
  double x1= 2*pow(t,2)/d;
  double x2= pow(t,2)/d;
  double y1=(8 - 2*m*pow(t,2))/d;
  double y2=(-4 - 4*rp*t - m*pow(t,2))/d;

  if(DEBUG) {
    printf("[call_205 #%d] t= %.10lf\n", P, t);
    printf("[call_205 #%d] w= %.10lf\n", P, w);
    printf("[call_205 #%d] m= %.10lf\n", P, m);
    printf("[call_205 #%d] d= %.10lf\n", P, d);
    printf("[call_205 #%d] x0= %.10lf\n", P, x0);
    printf("[call_205 #%d] x1= %.10lf\n", P, x1);
    printf("[call_205 #%d] x2= %.10lf\n", P, x2);
    printf("[call_205 #%d] y1= %.10lf\n", P, y1);
    printf("[call_205 #%d] y2= %.10lf\n", P, y2);
    printf("\n");
  }

  double k=0;
  if(LH==1) {
    k= -cos(w/2 + (double)1/2)/cos(w/2 - (double)1/2);
  }
  if(LH==0) {
    k=  sin((double)1/2 - w/2)/sin((double)1/2 + w/2);
  }
  if(DEBUG) printf("[call_205 #%d] k= %.10lf\n", P, k); 

  d= 1 + y1*k - y2*pow(k, 2);
  filter.a0= (x0 - x1*k + x2*pow(k,2))/d;
  filter.a1= (-2*x0*k + x1 +x1*pow(k,2) - 2*x2*k)/d;
  filter.a2= (x0*pow(k,2) -x1*k +x2)/d;
  filter.b1= (2*k + y1 + y1*pow(k,2) - 2*y2*k)/d;
  filter.b2= (-pow(k,2) - y1*k +y2)/d;

  if(DEBUG) {
    printf("[call_205 #%d] d= %.10lf\n", P, d);
    printf("[call_205 #%d] a0= %.10lf\n", P, filter.a0);
    printf("[call_205 #%d] a1= %.10lf\n", P, filter.a1);
    printf("[call_205 #%d] a2= %.10lf\n", P, filter.a2);
    printf("[call_205 #%d] b1= %.10lf\n", P, filter.b1);
    printf("[call_205 #%d] b2= %.10lf\n", P, filter.b2);
  }
  if(LH==1) {
    filter.a1=-filter.a1;
    filter.b1=-filter.b1;
  }

  if(DEBUG) {
    printf("[call_205 #%d] a1= %.10lf\n", P, filter.a1);
    printf("[call_205 #%d] b1= %.10lf\n", P, filter.b1);
  }

  return(filter);
}

/*
  NP : Number of poles / Order
  PR: Percent ripple
  LH: 0 for lowpass, 1 for highpass
  FC: cutoff frequency (0 to 0.5), typically fc(Hz)/SR(Hz), with fc<SR/2. Example: (double)1000/44100
*/
ChebFilter create_che_filter(int NP, double PR, int LH, double FC) {

  ChebFilter filter = {0, 0, 0, 0, 0, 0};

  double a[22];
  double b[22];

  double ta[22];
  double tb[22];

  for(int i=0; i<22; i++) {
    a[i]=0;
    b[i]=0;
  }

  a[2]=1;
  b[2]=1;


  for(int p=0; p<NP/2; p++) {
    filter=call_205(p, filter, FC, NP, LH, PR);

    for(int i=0; i<22; i++) {
      ta[i]=a[i];
      tb[i]=b[i];
    }

    for(int i=2; i<22; i++) {
      a[i]=filter.a0*ta[i] + filter.a1*ta[i-1] + filter.a2*ta[i-2];
      b[i]=   tb[i] - filter.b1*tb[i-1] - filter.b2*tb[i-2];
    }
  }

  b[2]=0;
  for(int i=0; i<20; i++) {
    a[i]=a[i+2];
    b[i]=-b[i+2];
  }

  double sa=0;
  double sb=0;

  for(int i=0; i<20; i++) {
    if(LH==0) { 
      sa=sa+a[i];
      sb=sb+b[i];
    }
    if(LH==1) { 
      sa=sa+a[i]*pow(-1, i);
      sb=sb+b[i]*pow(-1, i);
    }
  }

  double gain=sa/(1-sb);

  for(int i=0; i<20; i++) 
    a[i]=a[i]/gain;
  
  printf("Order: %d\n", NP);
  printf("Fc: %lf\n", FC);
  printf("a0: %lf\n", filter.a0);
  printf("a1: %lf\n", filter.a1);
  printf("a2: %lf\n", filter.a2);
  printf("b1: %lf\n", filter.b1);
  printf("b2: %lf\n", filter.b2);

  filter.x0=0;
  filter.x1=0;
  filter.x2=0;
  filter.y1=0;
  filter.y2=0;

  return(filter);
}

ChebFilter create_bw_low_pass_filter(int NP, double FC) {
  return(create_che_filter(NP, 0, 0, FC));
}

ChebFilter create_bw_high_pass_filter(int NP, double FC) {
  return(create_che_filter(NP, 0, 1, FC));
}

ChebFilter create_che_low_pass_filter(int NP, double FC, double PR) {
  return(create_che_filter(NP, PR, 0, FC));
}

ChebFilter create_che_high_pass_filter(int NP, double FC, double PR) {
  return(create_che_filter(NP, PR, 1, FC));
}

double bw_low_pass(ChebFilter filter, double x0) {
  double y=0;

  y= filter.a0*x0 + filter.a1*filter.x1 + filter.a2*filter.x2 + filter.b1*filter.y1 + filter.b2*filter.y2;

  filter.y2=filter.y1;
  filter.y1=y;
  filter.x2= filter.x1;
  filter.x1= filter.x0;
  filter.x0=x0;

  return(y);
}

int main() {
  ChebFilter filter = create_bw_low_pass_filter(4, (double)205/44100);

  printf("Filtering 0.5678: %.10lf\n", bw_low_pass(filter, 0.5678));
  return(0);
}