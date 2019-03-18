#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "includes/fraction.h"

long long pgcd(long long a, long long b)
{
	long long c;
    while(b!=0){
        c=a%b;
        a=b;
        b=c;
    }
    return a;
}

long long ppcm(long long a, long long b)
{
  return((a/pgcd(a,b)) * b);
}

fraction_t fraction_simplify(fraction_t fract)
{
	long long r=pgcd(fract.n, fract.d);
    fract.n/=r;
    fract.d/=r;
    if(fract.d<0){
        fract.n=-fract.n;
        fract.d=-fract.d;
	}
	return fract;
}

void fraction_swap(long long *a, long long *b)
{
    long long tmp = *a;
    *a = *b;
    *b = tmp;
}

fraction_t fraction_operator(fraction_t fract1, fraction_t fract2, char operator)
{
	long long pcm;
    fract1 = fraction_simplify(fract1);
    fract2 = fraction_simplify(fract2);
    switch(operator)
	{
		case '+': 
            pcm = ppcm(fract1.d, fract2.d);
            return(fraction_simplify((fraction_t){pcm/fract1.d*fract1.n + pcm/fract2.d*fract2.n, pcm}));
		case '-':
            pcm = ppcm(fract1.d, fract2.d);
            return(fraction_simplify((fraction_t){pcm/fract1.d*fract1.n - pcm/fract2.d*fract2.n, pcm}));	
		case '*': 
            fraction_swap(&fract1.n, &fract2.n);
            fract1 = fraction_simplify(fract1);
            fract2 = fraction_simplify(fract2);
            return(fraction_simplify((fraction_t){fract1.n * fract2.n, fract1.d * fract2.d}));	
		case '/': 
            fraction_swap(&fract2.n, &fract2.d);
            return(fraction_operator(fract1, fract2, '*'));	
	}
	return(FRACT_NULL);
}

int fraction_test(fraction_t fract1, fraction_t fract2, char *test)
{
    fract1 = fraction_simplify(fract1);
    fract2 = fraction_simplify(fract2);
    if(strcmp("==", test) == 0)return(fract1.n * fract2.d == fract2.n * fract1.d);
    if(strcmp(">", test) == 0)return(fract1.n * fract2.d > fract2.n * fract1.d);	
    if(strcmp(">=", test) == 0)return(fract1.n * fract2.d >= fract2.n * fract1.d);
    if(strcmp("<", test) == 0)return(fract1.n * fract2.d < fract2.n * fract1.d);
    if(strcmp("<=", test) == 0)return(fract1.n * fract2.d <= fract2.n * fract1.d);
	return 0;
}


fraction_t fraction_power(fraction_t fract, int lambda)
{
	long long i;
	fraction_t one = {1,1};
	if (lambda == 0) return one;
	if (lambda < 0) {
		i=fract.n;
		fract.n=fract.d;
		fract.d=i;
	}
	for(i=1; i<abs(lambda); i++){
		fract.n*=fract.n;
		fract.d*=fract.d;
	}
	return fraction_simplify(fract);
}

fraction_t double2fraction(double lambda){
	long long count = 1;
	double num = fabs(lambda);
	num -= (int)num;
	while(num != 0.0) {
		num *= 10;
		count *= 10;
		num -= (int)num;
	}
	fraction_t fraction = {(long long)(lambda*count), count};
	return fraction_simplify(fraction);
}

double str2double(const char *str, int len)
{
    char *s = calloc(len+1, sizeof(char));
    memcpy(s, str, len);
    double ret = strtod(s,NULL);
    free(s);
    return ret;
}

fraction_t string2fraction(const char *str, int len)
{
	int i;
	fraction_t num, den;
	for (i = 0; i < len; i++) {
		if (str[i] == '/') break;
	}
	num = double2fraction(str2double(str, i));
	if(i++ == len)return(num);
	den = double2fraction(str2double(str+i, len-i));
	return fraction_operator(num, den, '/');
}

char *fraction2string(fraction_t fract)
{
	char * s = malloc(50 * sizeof(char));
	if(fract.d == 1){
		sprintf(s, "%lld", fract.n);
	} else {
		sprintf(s, "%lld/%lld", fract.n, fract.d);
	}
	return s;
}