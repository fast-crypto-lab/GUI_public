

#include <stdio.h>

#include "gf2ext_poly.hpp"
#include "benchmark.h"

#include "stat_profile.h"

#include "config.h"

//#define EXT 127
//#define MAXD 9
//#define EXT 103
//#define MAXD 129
//#define EXT 128
//#define MAXD 9
#define EXT CORE_SIZE
#define MAXD MAX_DEG



#define TEST_RUN 1000



typedef GF2EXT<EXT> gf_t;
typedef poly<EXT,MAXD> poly_t;


unsigned deg( unsigned dd )
{
	if( 0 == dd ) return 0;
	if( 1 == dd ) return 1;
	if( 2 == dd ) return 1;
	if( 3 == dd ) return 2;
	if( 4 == dd ) return 1;

	unsigned tt = 1;

	while( dd >= tt ) tt <<= 1;
	tt >>= 1;

	return 1 + deg( dd-tt );
}

void rand_quad_poly( poly_t & p )
{
	gf_t r;
	p.set_zero();
	for(int i=MAXD;i>=0;i--) {
		if( 2 >= deg(i) ) {
			r = gf_t::rand();
			p.append_term( r , i );
		}
	}
}



int main()
{

	printf("====\ntest: gf(2^%d) poly: max deg: %d\n====\n",EXT,MAXD);
	srand(17);

	poly_t p;
	gf_t inp = gf_t::rand();
	rand_quad_poly( p );
	gf_t r = p.eval( inp );

	p += r;
	r = p.eval( inp );


	printf( "rand gf: " );
	inp.fdump( stdout );

	printf("\nrand poly: ");
	p.fdump2( stdout );
	printf("\neval: ");
	r.fdump( stdout );

	poly_t x_2ext_x;
	//x_2ext_x.append_term(gf_t::one(),2);
	//x_2ext_x.append_term(inp,1);
	p.X_2ext_X( x_2ext_x );
	printf("\n\nx_2ext_x: ");
	x_2ext_x.fdump2( stdout );

	poly_t gcd;
	poly_t::euclid_gcd( gcd , p , x_2ext_x );

	printf("\n\ngcd: ");
	gcd.fdump2( stdout );

	bool flag = p.find_unique_root( r );

	printf("\n\nfind unique root: %d , " , flag );
	r.fdump( stdout );

/////////////  benchmark //////////////////////

	printf("\n\nbenchamrk: .....%d\n", TEST_RUN );

	benchmark bm_ext;
	benchmark bm_gcd;
	bm_init( & bm_ext );
	bm_init( & bm_gcd );

#ifdef CONFIG_PROFILE
stat_reset();
#endif

	unsigned unique = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {

		rand_quad_poly( p );
		r = p.eval( inp );
		p += r;
#ifdef CONFIG_PROFILE
profile_start();
#endif
		BENCHMARK( bm_ext , {
		p.X_2ext_X( x_2ext_x );
		});

		BENCHMARK( bm_gcd , {
		poly_t::euclid_gcd( gcd , p , x_2ext_x );
		});
#ifdef CONFIG_PROFILE
profile_end();
#endif
		if( 1 == gcd.deg[0] ) unique++;
	}

	printf("unique root rate: %d/%d\n\n",unique,TEST_RUN);

	char msg[256];
	bm_dump( msg , 256 , &bm_ext );
	printf("poly extension:\n");
	puts(msg);

	bm_dump( msg , 256 , &bm_gcd );
	printf("\npoly gcd:\n");
	puts(msg);

	printf("\n\n");

#ifdef CONFIG_PROFILE
stat_say();
#endif

	return 0;
}
