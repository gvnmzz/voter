#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include"GSamp.h"

double smprng(void);

int main(void)
{
	int k, i, j, n=10, seq[10]={2,2,2,1,1,0,0,0,0,0};
	graph G;

	srand((long int) time(NULL));
	gsaminit(seq,n);

	printf("Sampling with uniform choice on the allowed nodes\n");
	for (k=0; k<5; k++) {
		G = gsam(smprng,0);
		printf("Adjacency list\n");
		for (i=0; i<n; i++) {
			printf("%d: ",i);
			for (j=0; j<seq[i]; j++) printf("%d ",G.list[i][j]);
			printf("\n");
		}
		printf("\n");
		printf("log(W): %.10f\n\n\n",G.weight);
	}

	printf("Sampling with uniform choice on the allowed stubs\n");
	for (k=0; k<5; k++) {
		G = gsam(smprng,1);
		printf("Adjacency list\n");
		for (i=0; i<n; i++) {
			printf("%d: ",i);
			for (j=0; j<seq[i]; j++) printf("%d ",G.list[i][j]);
			printf("\n");
		}
		printf("\n");
		printf("log(W): %.10f\n\n\n",G.weight);
	}

	gsamclean();

	return 0;
}

double smprng(void)
{
	return ((double) rand())/((double) RAND_MAX);
}
