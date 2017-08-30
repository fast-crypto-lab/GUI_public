
#include "run_config.h"

#include "gf2ext.hpp"
#include "benchmark.h"

#include <stdio.h>

#define TEST_NUM 10000



//#define TARGET "gf2ext_neon"
#define TARGET "gf2^128  sse"
//#define TARGET "gf2^94  sse"

//typedef gf2ext_u64<103> gf;
//typedef gf2ext_sse<103> gf;
//typedef gf2ext_sse<95> gf;
//typedef gf2ext_sse<94> gf;

void dump( const uint8_t *v )
{
	const uint32_t * v32 = (const uint32_t *)v;
	printf("%08x %08x, %08x %08x", v32[3] , v32[2] , v32[1] , v32[0] );
}


template <typename gf>
void bm( char * name )
{

	gf a[TEST_NUM],inv_a[TEST_NUM],check[TEST_NUM];

	char text[256];
	benchmark bm_mul,bm_inv,bm_squ;

	printf("================ benchmakr: %s  ===============\n",name);
	printf("sizeof(gf): %d\n", sizeof(gf));
	for (unsigned i=0; i<TEST_NUM; i++) {
		a[i] = gf::rand();
	}

	bm_init( & bm_inv );
	bm_start( & bm_inv );
	for(unsigned i=0;i<TEST_NUM;i++)
		inv_a[i] = a[i].inv();
	bm_stop( & bm_inv );
	bm_dump( text , 256 , &bm_inv );
	printf("benchmark inv():\n%s\n\n", text );

	bm_init( & bm_mul );
	bm_start( & bm_mul );
	for(unsigned i=0;i<TEST_NUM;i++)
		check[i] = a[i] * inv_a[i] ;
	bm_stop( & bm_mul );
	bm_dump( text , 256 , &bm_mul );
	printf("benchmark mul():\n%s\n\n", text );

	bm_init( & bm_squ );
	bm_start( & bm_squ );
	for(unsigned i=0;i<TEST_NUM;i++)
		check[i] = a[i].squ();
	bm_stop( & bm_squ );
	bm_dump( text , 256 , &bm_squ );
	printf("benchmark squ():\n%s\n\n", text );


}


int main()
{


	typedef gf2ext_u64<96> gf_t;
	gf_t a , b , c;
	a = gf_t::rand();
	b = a.inv();
	c = gf_t::rand();

	printf("a: "); a.fdump2( stdout );
	printf("\nb: "); b.fdump2( stdout );
	c = a*b;
	printf("\nc: "); c.fdump2( stdout );
	printf("\n\n");

	printf("Testing: %s\n",TARGET);
	printf("Run %d times.\n",TEST_NUM);

//typedef gf2ext_u64<103> gf;
//typedef gf2ext_sse<103> gf;
//typedef gf2ext_sse<95> gf;


//	bm< gf2ext_u64<94> >("gf2ext_u64<94>");
//	bm< gf2ext_u64<95> >("gf2ext_u64<95>");
//	bm< gf2ext_u64<96> >("gf2ext_u64<96>");
//	bm< gf2ext_u64<103> >("gf2ext_u64<103>");
	bm< gf2ext_u64<127> >("gf2ext_u64<127>");
	bm< gf2ext_u64<128> >("gf2ext_u64<128>");


#ifdef CONFIG_HAS_PCLMULQDQ
//	bm< gf2ext_sse<94> >("gf2ext_sse<94>");
//	bm< gf2ext_sse<95> >("gf2ext_sse<95>");
//	bm< gf2ext_sse<96> >("gf2ext_sse<96>");
//	bm< gf2ext_sse<103> >("gf2ext_sse<103>");
	bm< gf2ext_sse<127> >("gf2ext_sse<127>");
	bm< gf2ext_sse<128> >("gf2ext_sse<128>");
#endif

#ifdef CONFIG_PSHUFB
	bm< gf2ext_sse_tbl<96> >("gf2ext_sse_tbl<96>");
	bm< gf2ext_sse_tbl<128> >("gf2ext_sse_tbl<128>");
#endif

#ifdef CONFIG_NEON
	bm< gf2ext_neon<96> >("gf2ext_neon<96>");
	bm< gf2ext_neon<128> >("gf2ext_neon<128>");
#endif

	return 0;
}
