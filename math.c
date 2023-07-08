#include <stdio.h>
#include <math.h>

#define 	PR 10
#define		NP	4
#define 	P 2

int main() {

  double rp= -cos(M_PI/(NP*2) + (P-1)*M_PI/NP);
  double ip=  sin(M_PI/(NP*2) + (P-1)*M_PI/NP);

  double es= sqrt(pow(100.0 / (100.0-PR),2) - 1.0);

	// VX = (1/NP) * LOG( (1/ES) + SQR( (1/ES^2) + 1) )
	double vx= (1.0/NP) * log( (1.0/es) + sqrt( (1.0/pow(es,2))+1.0) );
	double kx= (1.0/NP) * log( (1.0/es) + sqrt( (1.0/pow(es,2))-1.0) );
				kx= (exp(kx) + exp(-kx))/2;
				rp= rp * ((exp(vx) - exp(-vx))/2.0)/kx;
				ip= ip * ((exp(vx) + exp(-vx))/2.0)/kx;

	printf("rp= %.10lf!\n", rp);
	printf("ip= %.10lf!\n", ip);
	printf("es= %.10lf!\n", es);
	printf("vx= %.10lf!\n", vx);
	printf("kx= %.10lf!\n", kx);

	return(0);
}