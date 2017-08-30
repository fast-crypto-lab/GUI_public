

#include "quartz_core.h"

#include "crand.h"
#include <stdio.h>

#include "benchmark.h"

#ifdef CONFIG_PROFILE
#include "stat_profile.h"
#endif



const unsigned ext = CORE_SIZE;
const unsigned max_d = MAX_DEG;
const unsigned minus = MINUS;
const unsigned vinegar = VINEGAR;

const unsigned num_n = ext+vinegar;
const unsigned num_m = ext-minus;


#define TEST_KEY 10
#define TEST_RUN 1000



int main()
{

	VEC<num_m> z0;
	VEC<num_m> z1;
	VEC<num_n> w0;
	VEC<num_n> w1;

	uint8_t sha_seed[_LEN_SHA256_] = {0};

	quartz_pub_key_t pub_key;
	quartz_sec_key_t sec_key;

	printf("benchmark (%d):  quartz<%d,%d,%d,%d>\n\n", TEST_RUN , ext , max_d , minus , vinegar );

	struct benchmark bm_sec;
	bm_init( & bm_sec );

	struct benchmark bm_pub;
	bm_init( & bm_pub );

	struct benchmark bm_key;
	bm_init( & bm_key );

	printf("benchmark genkey() (%d)\n", TEST_KEY);

	for(unsigned i=0;i<TEST_KEY;i++) {
BENCHMARK( bm_key , {
		quartz_gen_key(pub_key,sec_key);
} );
	}

#ifdef CONFIG_PROFILE
	stat_reset();
#endif

	printf("benchmark sign/verify (%d)\n", TEST_RUN);

	for(unsigned i=0;i<TEST_RUN;i++) {
		z0 = VEC<num_m>::rand();
		memset(sha_seed,0,_LEN_SHA256_);
		memcpy(sha_seed,&z0,(num_m+7)/8);
BENCHMARK( bm_sec , {
		quartz_sec_map( w0 , sec_key , z0 , sha_seed);
} );
BENCHMARK( bm_pub , {
		mpkc_pub_map( z1 , pub_key , w0 );
} );
	}

	printf("\n----------------------------\n\n");
	char msg[256];
	bm_dump( msg , 256 , & bm_key );
	printf("gen key():\n%s\n",msg);

	bm_dump( msg , 256 , & bm_sec );
	printf("sec map():\n%s\n",msg);

	bm_dump( msg , 256 , & bm_pub );
	printf("pub map():\n%s\n",msg);

#ifdef CONFIG_PROFILE
	bm_dump( msg , 256 , & prepare_poly );
	printf("prepare_poly():\n%s\n",msg);
	bm_dump( msg , 256 , & do_poly );
	printf("do_poly():\n%s\n",msg);
	bm_dump( msg , 256 , & ext_gcd_poly );
	printf("ext_gcd_poly():\n%s\n",msg);

//	stat_reset();

	printf("poly add(): [%ld,%ld,%ld]\n",num_poly(0),num_poly(1),num_poly(2));
	printf("gf mul(): [%ld,%ld,%ld]\n",num_mul(0),num_mul(1),num_mul(2));
	printf("gf squ(): [%ld,%ld,%ld]\n",num_squ(0),num_squ(1),num_squ(2));
	printf("gf inv(): [%ld,%ld,%ld]\n",num_inv(0),num_inv(1),num_inv(2));
#endif

	return 0;
}
