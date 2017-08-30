#include "stat_profile.h"

#include <stdint.h>

benchmark prepare_poly;
benchmark do_poly;
benchmark ext_gcd_poly;

unsigned __phase = 0;

#define NUM_PHASE 10

uint64_t __num_poly[NUM_PHASE] = {0};

uint64_t __num_mul[NUM_PHASE] = {0};

uint64_t __num_squ[NUM_PHASE] = {0};

uint64_t __num_inv[NUM_PHASE] = {0};


void stat_say()
{
	char text[256];

	printf("\n\n============ profile: ============\n");
	bm_dump( text , 256 , &prepare_poly );
	printf("prepare_poly():\n%s\n", text );
	bm_dump( text , 256 , &do_poly );
	printf("do_poly():\n%s\n", text );
	bm_dump( text , 256 , &ext_gcd_poly );
	printf("ext_gcd_poly():\n%s\n", text );

	printf("---\n");
	for(unsigned i=0;i<NUM_PHASE;i++) {
		printf("mul[%d]: %ld\n",i,num_mul(i));
		printf("squ[%d]: %ld\n",i,num_squ(i));
		printf("inv[%d]: %ld\n",i,num_inv(i));
		printf("poly[%d]: %ld\n",i,num_poly(i));
		printf("---\n");
	}
}

void profile_start()
{
	__phase = 1;
}

void profile_end()
{
	__phase = 0;
}


void new_phase()
{
	__phase++;
	if( NUM_PHASE == __phase ) __phase = 0;
}

uint64_t num_mul(unsigned phase)
{
	return __num_mul[phase+1];
}

uint64_t num_squ(unsigned phase)
{
	return __num_squ[phase+1];
}

uint64_t num_inv(unsigned phase)
{
	return __num_inv[phase+1];
}

uint64_t num_poly(unsigned phase)
{
	return __num_poly[phase+1];
}

void stat_reset()
{
	bm_init(&prepare_poly);
	bm_init(&do_poly);
	bm_init(&ext_gcd_poly);

	__phase = 0;
	for(int i=0;i<NUM_PHASE;i++) __num_mul[i]=0;
	for(int i=0;i<NUM_PHASE;i++) __num_squ[i]=0;
	for(int i=0;i<NUM_PHASE;i++) __num_inv[i]=0;
	for(int i=0;i<NUM_PHASE;i++) __num_poly[i]=0;
}

void count_mul()
{
	__num_mul[__phase]++;
}

void count_squ()
{
	__num_squ[__phase]++;
}

void count_inv()
{
	__num_inv[__phase]++;
}

void count_poly_add()
{
	__num_poly[__phase] ++;
}



