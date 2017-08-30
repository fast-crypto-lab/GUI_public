#include "blas.h"

#include <stdio.h>
#include <stdlib.h>


typedef mat_u64<127,127> mat;
typedef VEC<127> vec;

int main()
{
	mat a;
	mat b;
	vec vv;
	vec vv2;

	uint8_t msg[16] = {0};

	bool good = mat::rand_inv( a , b );
	printf("rand(): %d\n",good);

	mat c;
	mat::mul( c , a , b );

	printf("c:\n");
	c.fdump(stdout);

	printf("\n\n");

	vv = vec::rand();
	printf("vv: "); vv.fdump(stdout); printf("\n");

	vv2 = a.prod( vv );
	printf("[a] x vv: "); vv.fdump(stdout); printf("\n");

	vv2 = b.prod( vv2 );
	printf("[b] x ?1: "); vv.fdump(stdout); printf("\n");

	vv2 ^= vv;
	printf("test pass? %d\n", vv2.is_zero() );

	return 0;

}


