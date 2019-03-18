#ifndef FRACTION
#define FRACTION

#define FRACT_UNDEFINED (fraction_t){-1,1}
#define FRACT_NULL (fraction_t){0,1}
#define FRACT_ONE (fraction_t){1,1}
typedef struct {
    long long n;
    long long d;
} fraction_t;

long long pgcd(long long a, long long b);
long long ppcm(long long a, long long b);
void fraction_swap(long long *a, long long *b);
fraction_t fraction_simplify(fraction_t fract);
fraction_t fraction_operator(fraction_t fract1, fraction_t fract2, char operator);
int fraction_test(fraction_t fract1, fraction_t fract2, char *test);
fraction_t fraction_power(fraction_t fract, int lambda);
fraction_t double2fraction(double lambda);
fraction_t string2fraction(const char *str, int len);
char * fraction2string(fraction_t fract);

#endif