#include <stdio.h>
#include <stdlib.h>
#define MATHLIB_STANDALONE 1
#include <Rmath.h>


int main(int argc, char **argv)
{
		
	double x = atof(argv[1]);
	double mu = atof(argv[2]);
	double lambda = atof(argv[3]);
	
	double leftarg = sqrt(lambda/x)*(x/mu - 1.0);
	double middlearg = (2*lambda)/mu;
	double rightarg = -(sqrt(lambda/x)*(x/mu + 1.0));
	
	double newmu = 0;
	double newsigma = 1;
	int lower_tail =1;
	int giveLog = 1;
	double left = pnorm(leftarg, newmu, newsigma, lower_tail, giveLog);
	
	giveLog = 1;
	double right = pnorm(rightarg, newmu, newsigma, lower_tail, giveLog);
	
	double pval = left + exp(middlearg + right);
	//printf("%f\n", pval);
	printf("Leftarg: %f, Middlearg: %f, Rightarg: %f\n", left, middlearg, right);
	return pval;
}
