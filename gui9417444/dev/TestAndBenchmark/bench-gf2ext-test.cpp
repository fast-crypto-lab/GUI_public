

#include <stdio.h>

#include "benchmark.h"


#define TEST_NUM 10000

#include "gf2ext.hpp"
#include "gf2ext-sse.hpp"




template <typename gf_t>
void test_gf( const char * name )
{
	printf("\n==============  testing: %s  ============\n",name);
	printf("sizeof(gf_t): %u\n", sizeof(gf_t));

	gf_t a,inv_a,check;

	printf("a: "); a.fdump( stdout ); printf("\n");
	printf("inv_a: "); inv_a.fdump( stdout ); printf("\n");
	printf("check: "); check.fdump( stdout ); printf(" : is_one(): %d , is_zero(): %d\n", check.is_one() , check.is_zero() );
	a = gf_t::rand();
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

	struct benchmark bmul,bsqu,binv;
	bm_init( & bmul );
	bm_init( & bsqu );
	bm_init( & binv );
	BENCHMARK( bmul , {
		for(unsigned i=TEST_NUM;i>0;i--) a *= a;
	} );
	BENCHMARK( bsqu , {
		for(unsigned i=TEST_NUM;i>0;i--) a = a.squ();
	} );
	BENCHMARK( binv , {
		for(unsigned i=TEST_NUM;i>0;i--) a = a.inv();
	} );

	char msg[256];
        bm_dump( msg , 256 , & bmul );
        printf("mul():\n%s\n",msg);
        bm_dump( msg , 256 , & bsqu );
        printf("squ():\n%s\n",msg);
        bm_dump( msg , 256 , & binv );
        printf("inv():\n%s\n",msg);

}



int main()
{

	test_gf<gf2ext_u64<94> >( "gf2ext_u64<94>" );
	test_gf<gf2ext_u64<95> >( "gf2ext_u64<95>" );
	test_gf<gf2ext_u64<96> >( "gf2ext_u64<96>" );
	test_gf<gf2ext_u64<103> >( "gf2ext_u64<103>" );
	test_gf<gf2ext_u64<127> >( "gf2ext_u64<127>" );

	test_gf<gf2ext_sse<94> >( "gf2ext_sse<94>" );
	test_gf<gf2ext_sse<95> >( "gf2ext_sse<95>" );
	test_gf<gf2ext_sse<96> >( "gf2ext_sse<96>" );
	test_gf<gf2ext_sse<103> >( "gf2ext_sse<103>" );
	test_gf<gf2ext_sse<127> >( "gf2ext_sse<127>" );

/*
	a = gf::one();
	//a.v = _mm_slli_epi32(a.v,1);
	a.v = _mm_slli_si128(a.v,8);
	inv_a = a;
	check = a * inv_a;
	printf("a   : "); a.fdump2( stdout ); printf("\n");
	printf("b   : "); inv_a.fdump2( stdout ); printf("\n");
	printf("axb : "); check.fdump2( stdout ); printf("\n");

	a = gf::irrPoly();
	a.v = _mm_slli_epi64(a.v,34);
	printf("a   : "); a.fdump2( stdout ); printf("\n");

	a = gf::one();
	a.v = _mm_slli_epi32(a.v,1);
	inv_a = gf::rand();
	check = a * inv_a;
	printf("a   : "); a.fdump2( stdout ); printf("\n");
	printf("b   : "); inv_a.fdump2( stdout ); printf("\n");
	printf("axb : "); check.fdump2( stdout ); printf("\n");

	a = gf::one();
	//a.bit_shl(1);
	inv_a = a.inv();
	check = a * inv_a;
*/

	return 0;
}
