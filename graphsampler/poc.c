#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "PLSG.h"
#include "GSamp.h"
#include "mt64.h"


double smprng(void);

int main(void)
{
	int i, j, k, kn, n=1000, *seq;
	int nswp = 10000;
	graph G;
	double op[n], dist;

	FILE *f = fopen("file.txt", "w");

	init_genrand64(time(NULL));

	srand((long int) time(NULL));

	seq = plseqgen(n, 3.5, smprng);
	printf("First sequence. 1000 nodes, gamma = 2.5\n");
	for (i=0; i<n-1; i++) printf("%d, ",seq[i]);
	printf("%d\n\n",seq[99]);

	gsaminit(seq,n);

	printf("Sampling with uniform choice on the allowed nodes\n");
	for (k=0; k<1; k++) {
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

// Now define the vector of opinions
for (i=0; i<n-1;i++) {
	op[i] = genrand64_real2();
}

// Now we want to make it evolve
for (i=0; i<nswp;i++) {
	for (j=0;j<n;j++) {
		// Pick a random node
		k = floor(genrand64_real2()*n);
		// Pick a random neighbour
		kn = ceil(genrand64_real2()*seq[k]);
		//dist = fabs(op[k]-op[G.list[k][kn]]);
    //if (dist < 0.2) {
			op[k] = 0.8*op[k] + 0.2*op[G.list[k][kn]];
			fprintf(f,"%d %f \n",i,op[k]);
		//}
    //printf("%d %d %d %f\n",j,kn,G.list[k][kn],dist);
	}
}





	free(seq);
	return EXIT_SUCCESS;
}

double smprng(void)
{
	return ((double) rand())/((double) RAND_MAX);
}
