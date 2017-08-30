

#include <stdio.h>

#define TEST_NUM 1000

#include "gf2ext.hpp"

#include "run_config.h"

#define TEST_U64


template <typename gf_t>
void test_gf( const char * name )
{
	printf("\n==============  testing: %s  ============\n",name);
	printf("sizeof(gf_t): %u\n", sizeof(gf_t));

	gf_t a,inv_a,check;

	printf("a: "); a.fdump( stdout ); printf("\n");
	printf("inv_a: "); inv_a.fdump( stdout ); printf("\n");
	printf("check: "); check.fdump( stdout ); printf(" : is_one(): %d , is_zero(): %d\n", check.is_one() , check.is_zero() );

	a = gf_t::one();
	inv_a = a.inv();
	check = a * inv_a;
	printf("a: "); a.fdump( stdout ); printf("\n");
	printf("inv_a: "); inv_a.fdump( stdout ); printf("\n");
	printf("check: "); check.fdump( stdout ); printf(" : is_one(): %d\n", check.is_one() );

	for (unsigned i=0; i<TEST_NUM; i++) {
		a = gf_t::rand();
		inv_a = a.inv();
		check = a * inv_a;

		if( ! check.is_one() ) {

	printf("a: "); a.fdump( stdout ); printf("\n");
	printf("inv_a: "); inv_a.fdump( stdout ); printf("\n");
	printf("check: "); check.fdump( stdout ); printf(" : is_one(): %d\n", check.is_one() );

			printf("test fail (%d)!!!\n",i+1);
			return;
		}
	}

	printf("test passed(%d)!\n", TEST_NUM);

}



int main()
{

	uint8_t aa[16] = {0x10,0x32,0x54,0x76, 0x98,0xba,0xdc,0xfe, 0xf8,0xf9,0xf1,0xcc, 0,0,0,0};
	uint8_t bb[16] = {0x01,2,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

	uint32_t cc[4] = {0};

/*
	typedef gf2ext_neon<96> gf_t;

	gf_t a(aa);
	a.fdump2(stdout);
	printf("\n");

	gf_t b(bb);
	b.fdump2(stdout);
	printf("\n");

	gf_t c = a*b;
	c.fdump2(stdout);
	printf("\n");

	gf_t d = b.squ();
	d.fdump2(stdout);
	printf("\n");
*/


#ifdef TEST_U64
	test_gf<gf2ext_u64<94> >( "gf2ext_u64<94>" );
	test_gf<gf2ext_u64<95> >( "gf2ext_u64<95>" );
	test_gf<gf2ext_u64<96> >( "gf2ext_u64<96>" );
	test_gf<gf2ext_u64<103> >( "gf2ext_u64<103>" );
	test_gf<gf2ext_u64<127> >( "gf2ext_u64<127>" );
	test_gf<gf2ext_u64<128> >( "gf2ext_u64<128>" );
#endif

#ifdef CONFIG_HAS_PCLMULQDQ
	test_gf<gf2ext_sse<94> >( "gf2ext_sse<94>" );
	test_gf<gf2ext_sse<95> >( "gf2ext_sse<95>" );
	test_gf<gf2ext_sse<96> >( "gf2ext_sse<96>" );
	test_gf<gf2ext_sse<103> >( "gf2ext_sse<103>" );
	test_gf<gf2ext_sse<127> >( "gf2ext_sse<127>" );
	test_gf<gf2ext_sse<128> >( "gf2ext_sse<128>" );
#endif

#ifdef CONFIG_PSHUFB
	test_gf<gf2ext_sse_tbl<96> >( "gf2ext_sse_tbl<96>" );
	test_gf<gf2ext_sse_tbl<128> >( "gf2ext_sse_tbl<128>" );
#endif

#ifdef CONFIG_NEON
	test_gf<gf2ext_neon<96> >( "gf2ext_neon<96>" );
	test_gf<gf2ext_neon<128> >( "gf2ext_neon<128>" );
#endif
	return 0;
}
