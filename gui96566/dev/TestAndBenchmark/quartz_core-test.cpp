

#include "quartz_core.h"
#include <string.h>


const unsigned ext = CORE_SIZE;
const unsigned max_d = MAX_DEG;
const unsigned minus = MINUS;
const unsigned vinegar = VINEGAR;

const unsigned num_n = ext+vinegar;
const unsigned num_m = ext-minus;

#define TEST_RUN 100



int main()
{

	uint8_t sha_seed[_LEN_SHA256_] = {0};

	VEC<num_m> z0;
	VEC<num_m> z1;
	VEC<num_n> w0;
	VEC<num_n> w1;

	srand(17);

	quartz_pub_key_t pub_key;
	quartz_sec_key_t sec_key;

	quartz_gen_key(pub_key,sec_key);

	printf("test:  quartz<%d,%d,%d,%d>\n\n", ext , max_d , minus , vinegar );
	printf("size of pub_key: %d\n",pub_key.num_byte());
	printf("size of sec_key: %d\n",sec_key.num_byte());

	bool checked = true;

	z0 = VEC<num_m>::rand();
	printf("message: z0:");
	z0.fdump(stdout);
	printf("\n");

	memset(sha_seed,0,_LEN_SHA256_);
	memcpy(sha_seed,&z0,(num_m+7)/8);

	quartz_sec_map( w0 , sec_key , z0 , sha_seed);
	printf("signature: w0:");
	w0.fdump(stdout);
	printf("\n");

	mpkc_pub_map( z1 , pub_key , w0 );
	printf("check: z1:");
	z1.fdump(stdout);
	printf("\n");

	checked = (z0 == z1);
	printf("passed? %d\n", checked);

	printf("\nRunning %d tests: ", TEST_RUN);

	checked = true;
	for(unsigned i=0;i<TEST_RUN;i++) {

		z0 = VEC<num_m>::rand();

		memset(sha_seed,0,_LEN_SHA256_);
		memcpy(sha_seed,&z0,(num_m+7)/8);
		quartz_sec_map( w0 , sec_key , z0 ,sha_seed );
		mpkc_pub_map( z1 , pub_key , w0 );

		if( !(z0==z1) ) checked=false;
	}

	printf("passed %d tests? %d\n", TEST_RUN, checked);

	printf("\n\n");


	return 0;
}
