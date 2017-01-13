/*
   A C-program for MT19937
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)


   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdio.h>

/* Parametri */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* vectorul a */
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL

static unsigned long x[N]; /* tabloul de N numere */
static int mti=N+1; /* mti=N+1 inseamna ca  x[N] nu este initializat */

/* initializarea celor N numere ce se stocheaza in x[N] */

void init_genrand(unsigned long s)
{  /* modul de initializare a fost sugerat de Donald Knuth*/
    x[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++)
     {
        x[mti] =
	    (1812433253UL * (x[mti-1] ^ (x[mti-1] >> 30)) + mti);

        x[mti] &= 0xffffffffUL;

    }
}


/* Genereaza un numar aleator intreg din intervalul [0,0xffffffff] */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long v[2]={0x0UL, MATRIX_A};


    if (mti >= N)
       { /* generate N numere simultan */
        int j;

	if (mti == N+1)   /* daca init_genrand() nu a fost apelat, */
	    init_genrand(5489UL); /* un seed implicit */

        for (j=0;j<N-M;j++)
            {
            y = (x[j]&UPPER_MASK)|(x[j+1]&LOWER_MASK);
            x[j] = x[j+M] ^ (y >> 1) ^ v[y & 0x1UL];
            }
        for (;j<N-1;j++)
        {
            y = (x[j]&UPPER_MASK)|(x[j+1]&LOWER_MASK);
            x[j] = x[j+(M-N)] ^ (y >> 1) ^ x[y & 0x1UL];
        }
        y = (x[N-1]&UPPER_MASK)|(x[0]&LOWER_MASK);
        x[N-1] = x[M-1] ^ (y >> 1) ^ v[y & 0x1UL];

        mti = 0;
    }

    y = x[mti++];

    /* Operatiile de temperare */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}


/* generator de numere din intervalul inchis [0,1] */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* 4294967295= 2^32-1 */
}

/* generator de numere din intervalul [0,1) */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* 4294967296= 2^32 */
}

/* generator de numere din intervalul (0,1) */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);

}



