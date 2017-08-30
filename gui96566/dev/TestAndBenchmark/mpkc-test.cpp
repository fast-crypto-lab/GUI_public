
#include <stdio.h>

#include "mpkc.hpp"

#include "benchmark.h"

#define TEST_RUN 10000

#define N 95
#define M 100

const unsigned qterms = TERMS_QUAD_POLY(N);

typedef MAT<qterms,M> qpoly_t;

typedef VEC<M> vec_out_t;
typedef VEC<N> vec_in_t;


void quad_poly_eval(void * r ,const void * key,const void * i ){
	mpkc_pub_map<N,M>( *(VEC<M>*)r , *(const qpoly_t *)key , *(const VEC<N> *)i );
}

int main()
{
	qpoly_t poly;

	for(unsigned i=0;i<qterms;i++) poly.v[i] = vec_out_t::rand();

	qpoly_t poly2;
	interpolate<N,M>( poly2 , quad_poly_eval , &poly );

	bool checked = true;
	for(unsigned i=0;i<qterms;i++) {
		vec_out_t a = poly.v[i];
		a ^= poly2.v[i];

		if( ! a.is_zero() ) {
			printf("poly.v[%d] fail!\n",i);
			checked = false;
			break;
		}
	}

	printf("interpolate check: %s\n\n", checked?"success!":"fail!") ;



	vec_in_t inp;
	vec_out_t out;
	benchmark bm_eval;
	bm_init( & bm_eval );

	for(unsigned i=0;i<TEST_RUN;i++) {
		vec_out_t a;
		inp = vec_in_t::rand();
		bm_start( &bm_eval );
		mpkc_pub_map<N,M>( a , poly2 , inp );
		bm_stop( &bm_eval );
		out ^= a;
	}
	out.fdump(stdout);

	printf("\n\nbenchamrk: public poly evaluation:\n\n");
	char msg[256];
	bm_dump( msg , 256 , &bm_eval );

	printf("%s\n\n",msg);


	return 0;
}
