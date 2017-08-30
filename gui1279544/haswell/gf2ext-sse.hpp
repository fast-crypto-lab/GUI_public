#ifndef _GF2EXT_SSE_HPP_
#define _GF2EXT_SSE_HPP_

#include <stdint.h>



#ifdef CONFIG_PROFILE
#include "stat_profile.h"
#endif

#include "gf2ext.hpp"



#if defined( CONFIG_HAS_PCLMULQDQ )||defined( CONFIG_PSHUFB )
#include "emmintrin.h"
#endif

#if defined( CONFIG_PSHUFB )
#include "tmmintrin.h"
#endif

#if defined( CONFIG_HAS_PCLMULQDQ )
#include "wmmintrin.h"
#endif

#if defined( CONFIG_HAS_PCLMULQDQ )


template <unsigned W> struct gf2ext_sse;

template <unsigned W> inline
gf2ext_sse<W> _mul( const gf2ext_sse<W> & a , const gf2ext_sse<W> & b );



template <unsigned W>
struct gf2ext_sse {
	__m128i v;

	typedef gf2ext_sse gf_t;
	typedef gf2ext_u64<W> gf_u64;

	gf2ext_sse() { v^=v; }
	explicit gf2ext_sse(bool b) { *((gf_u64*)this)=gf_u64(b); } //{ int64_t b0=b; __m128i v0=_mm_cvtsi64_si128(b0); v0=_mm_shuffle_epi32(v0,0); v =_mm_cmpgt_epi32(v0,v0^v0); }
	gf2ext_sse( __m128i a ): v(a) {}
	gf2ext_sse( const gf_t & a ): v(a.v) {}
	gf2ext_sse( const gf_u64 & a ) { *((gf_u64*)this)=a; }
	explicit gf2ext_sse( const uint8_t * x ) { *((gf_u64*)this)=gf_u64(x); }

	gf_t squ() const;
	gf_t inv() const { return fermat_inv<gf_t,W>::exec(*this); }

	gf_t & operator *= ( const gf_t & a ) { *this = _mul( *this, a ); return *this; }
	gf_t & operator ^= ( const gf_t & a ) { v^=a.v; return *this; }
	gf_t & operator &= ( const gf_t & a ) { v&=a.v; return *this; }

	gf_t operator *( const gf_t & a ) const { return _mul(*this,a); }
	gf_t operator ^( const gf_t & a ) const { gf_t r=*this; r^=a; return r; }
	gf_t operator &( const gf_t & a ) const { gf_t r=*this; r&=a; return r; }

	gf_t & set_zero() { v^=v; return *this; }
	bool is_zero() const { return 0xffff==_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_setzero_si128(),v)); }
	//{ return ((const gf_u64 *)this)->is_zero(); }
	bool is_one() const { return ((const gf_u64 *)this)->is_one(); }

	static gf_t rand() { gf_t r = gf_u64::rand(); return r; }

	static const gf_t & zero() { return *(const gf_t *)&gf_u64::zero(); }
	static const gf_t & one() { return *(const gf_t *)&gf_u64::one(); }
	static const gf_t & irrPoly() { return *(const gf_t *)&gf_u64::irrPoly(); }

	static gf_t assign( const uint8_t * x ) { return gf_t(x); }
	static unsigned num_byte() { return gf_u64::num_byte(); }
	void dump( uint8_t * x ) const { ((const gf_u64 *)this)->dump(x); }
	void fdump(FILE *fp) const { ((const gf_u64 *)this)->fdump(fp); }
	void fdump2(FILE *fp) const { ((const gf_u64 *)this)->fdump2(fp); }

//	static const uint32_t _mask_56bit[4];

	static const uint32_t _mask_Wbit[4];
	static const uint64_t _reducer_W[2];


};


//__m128i _mm_clmulepi64_si128(__m128i v1, __m128i v2, const int imm8);


////////////////////////////////////////////////////////

template <unsigned W> inline 
void _reduce( __m128i &a0b0 , __m128i a1b1 );



template <unsigned W> inline
gf2ext_sse<W> gf2ext_sse<W>::squ() const
{
#ifdef CONFIG_PROFILE
	count_squ();
#endif
	__m128i a0b0 = _mm_clmulepi64_si128( v , v , 0x00 );
	__m128i a1b1 = _mm_clmulepi64_si128( v , v , 0x11 );

	_reduce<W>( a0b0 , a1b1 );

	gf_t r;
	r.v = a0b0;
	return r;
}



template <unsigned W> inline
gf2ext_sse<W> _mul( const gf2ext_sse<W> & a, const gf2ext_sse<W> & b )
{
#ifdef CONFIG_PROFILE
	count_mul();
#endif
	__m128i a0b0 = _mm_clmulepi64_si128( a.v , b.v , 0x00 );
	__m128i a1b1 = _mm_clmulepi64_si128( a.v , b.v , 0x11 );
// cross terms
#ifdef CONFIG_FAST_PCLMULQDQ
	__m128i a0b1 = _mm_clmulepi64_si128( a.v , b.v , 0x10 ) ^ _mm_clmulepi64_si128( a.v , b.v , 0x01 );
#else
	__m128i a0b1 = (__m128i)_mm_shuffle_pd( (__m128d)a.v,(__m128d)b.v,1);  /// --> [b_low64,a_high64]
	a0b1 = _mm_clmulepi64_si128( a0b1^a.v,a0b1^b.v, 0x10 ) ^ a0b0 ^ a1b1;
#endif
	a0b0 ^= (__m128i)_mm_shuffle_pd( (__m128d)_mm_setzero_si128() , (__m128d)a0b1 , 0 ); /// --> [a0b1_low,0]
	a1b1 ^= (__m128i)_mm_shuffle_pd( (__m128d)a0b1 , (__m128d)_mm_setzero_si128() , 1 ); /// --> [0,a0b0_high]
// end cross terms

	_reduce<W>( a0b0 , a1b1 );

	gf2ext_sse<W> r;
	r.v = a0b0;
	return r;
}




////////////////////////////    gf2103   reduce   ////////////////////////////////////////



inline void reduce_clmul_103( __m128i &a0b0 , __m128i a1b1 )  /// test fail. check!!!
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_mask_Wbit );
	__m128i rder = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_reducer_W );

	a0b0 ^= _mm_clmulepi64_si128( a1b1 , rder , 0x00 )
		^ _mm_slli_si128( _mm_clmulepi64_si128( a1b1 , rder , 0x01) , 8);
	a0b0 ^= _mm_clmulepi64_si128( _mm_srli_epi64(a0b0,39) , rder , 0x11 );
	//a0b0 ^= _mm_clmulepi64_si128( _mm_andnot_si128(msk,a0b0) , rder , 0x11 );
	a0b0 &= msk;

}


inline void reduce_shift( __m128i &a0b0 , __m128i a1b1 )
{
	//__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_mask_56bit );
	// {0xffffffff,0x00ffffff,0,0};
	__m128i msk = _mm_set_epi32( 0 , 0 , 0x00ffffff , 0xffffffff );
	__m128i a1b1h = _mm_andnot_si128( msk , a1b1 );
	a1b1 &= msk;

	__m128i shr103 = _mm_slli_epi64(a1b1,1);
	__m128i shr94 = _mm_slli_epi64(a1b1,2);
	shr103 = _mm_slli_si128(shr103,3);
	shr94 = _mm_slli_si128(shr94,4);
	a0b0 ^= shr103;
	a0b0 ^= shr94;
	shr103 = _mm_slli_si128(a1b1h,3);
	shr94 = _mm_slli_si128(a1b1h,4);
	shr103 = _mm_slli_epi64(shr103,1);
	shr94 = _mm_slli_epi64(shr94,2);
	a0b0 ^= shr103;
	a0b0 ^= shr94;

	msk = _mm_load_si128( (__m128i*) gf2ext_sse<103>::_mask_Wbit );
	__m128i a0b0h = _mm_andnot_si128( msk , a0b0 );
	a0b0 &= msk;
	a0b0h = _mm_srli_si128(a0b0h,12);
	a0b0 ^= _mm_srli_epi64(a0b0h,7);
	a0b0 ^= _mm_slli_epi64(a0b0h,2);
}


template <> inline 
void _reduce<103>( __m128i &a0b0 , __m128i a1b1 )
{
#ifdef CONFIG_PCLMULQDQ_REDUCE
	reduce_clmul_103( a0b0 , a1b1 );
#else
	reduce_shift( a0b0 , a1b1 );
#endif
}


////////////////////////////    gf296   reduce   ////////////////////////////////////////


//extern const uint32_t _mask_96bit[];
//extern const uint64_t _reducer_96[];

///  x^96 + x^10 + x^9 + x^6 + 1
///  x^128 + x^42 + x^41 + x^38 + x^32

/// [----- 192 ---- 128 ][------ 96 ----------]


inline void _reduce_96_pclmul( __m128i &a0b0 , __m128i a1b1 )
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<96>::_mask_Wbit );
	__m128i rder = _mm_load_si128( (__m128i*) gf2ext_sse<96>::_reducer_W );

	a0b0 ^= _mm_clmulepi64_si128( a1b1 , rder , 0x00 );
	__m128i a0b0h = _mm_srli_si128( a0b0 , 4 );
	a0b0 ^= _mm_clmulepi64_si128( a0b0h , rder , 0x11 );
	a0b0 &= msk;
}


inline void _reduce_96_sse( __m128i &a0b0 , __m128i a1b1 )
{
	_reduce_96_pclmul( a0b0 , a1b1 );
#if 0
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<95>::_mask_Wbit );
	__m128i msk32bit = _mm_set_epi32(0,0,0,0xffffffff);

	__m128i _x33 = _mm_slli_si128( _mm_slli_epi64(a1b1,1) , 4 );
	__m128i _x44low = _mm_slli_si128( _mm_slli_epi64(a1b1&msk32bit,4) , 5 );
	__m128i _x44high = _mm_slli_epi64( _mm_slli_si128( _mm_andnot_si128(msk32bit,a1b1) ,5), 4 );

	a0b0 ^= _x44high;
	__m128i _a0b0h = _mm_andnot_si128( msk , a0b0 );
	__m128i _x11 = _mm_srli_si128( _mm_srli_epi64(_a0b0h,4),10);
	__m128i _x0 = _mm_srli_si128( _mm_srli_epi64(_a0b0h,7),11);

	a0b0 ^= _x33 ^ _x44low ^ _x11 ^ _x0;
	a0b0 &= msk;
#endif
}



template <> inline
void _reduce<96>( __m128i &a0b0 , __m128i a1b1 )
{
//#ifdef CONFIG_FAST_PCLMULQDQ
#ifdef CONFIG_PCLMULQDQ_REDUCE
	_reduce_96_pclmul( a0b0 , a1b1 );
#else
	_reduce_96_sse( a0b0 , a1b1 );
#endif
}




////////////////////////////    gf295   reduce   ////////////////////////////////////////


//extern const uint32_t _mask_95bit[];
//extern const uint64_t _reducer_95[];


/// x^95 + x^11 + 1
/// x^128 + x^44 + x^33

inline void _reduce_95_sse( __m128i &a0b0 , __m128i a1b1 )
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<95>::_mask_Wbit );
	__m128i msk32bit = _mm_set_epi32(0,0,0,0xffffffff);

	__m128i _x33 = _mm_slli_si128( _mm_slli_epi64(a1b1,1) , 4 );
	__m128i _x44low = _mm_slli_si128( _mm_slli_epi64(a1b1&msk32bit,4) , 5 );
	__m128i _x44high = _mm_slli_epi64( _mm_slli_si128( _mm_andnot_si128(msk32bit,a1b1) ,5), 4 );

	a0b0 ^= _x44high;
	__m128i _a0b0h = _mm_andnot_si128( msk , a0b0 );
	__m128i _x11 = _mm_srli_si128( _mm_srli_epi64(_a0b0h,4),10);
	__m128i _x0 = _mm_srli_si128( _mm_srli_epi64(_a0b0h,7),11);

	a0b0 ^= _x33 ^ _x44low ^ _x11 ^ _x0;
	a0b0 &= msk;
}

inline void _reduce_95_pclmul( __m128i &a0b0 , __m128i a1b1 )
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<95>::_mask_Wbit );
	__m128i rder = _mm_load_si128( (__m128i*) gf2ext_sse<95>::_reducer_W );

	a0b0 ^= _mm_clmulepi64_si128( a1b1 , rder , 0x00 );
	__m128i a0b0h = _mm_srli_epi64( a0b0 , 31 );
	a0b0 ^= _mm_clmulepi64_si128( a0b0h , rder , 0x11 );
	a0b0 &= msk;
}


template <> inline
void _reduce<95>( __m128i &a0b0 , __m128i a1b1 )
{
//#ifdef CONFIG_FAST_PCLMULQDQ
#ifdef CONFIG_PCLMULQDQ_REDUCE
	_reduce_95_pclmul( a0b0 , a1b1 );
#else
	_reduce_95_sse( a0b0 , a1b1 );
#endif
}




////////////////////////////    gf294   reduce   ////////////////////////////////////////


//extern const uint32_t _mask_94bit[];
//extern const uint64_t _reducer_94[];

/// x^94 + x^21 + 1
/// x^128 + x^55 + x^34

inline void _reduce_94_sse( __m128i &a0b0 , __m128i a1b1 )
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<94>::_mask_Wbit );
	__m128i msk32bit = _mm_set_epi32(0,0,0,0xffffffff);

	__m128i _x34 = _mm_slli_si128( _mm_slli_epi64(a1b1,2) , 4 );
	__m128i _x55low = _mm_slli_si128( _mm_slli_epi64(a1b1&msk32bit,7) , 6 );
	__m128i _x55high = _mm_slli_epi64( _mm_slli_si128( _mm_andnot_si128(msk32bit,a1b1) ,6), 7 );

	a0b0 ^= _x55high;
	__m128i _a0b0h = _mm_andnot_si128( msk , a0b0 );
	__m128i _x21 = _mm_srli_si128( _mm_srli_epi64(_a0b0h,1),9);
	__m128i _x0 = _mm_srli_si128( _mm_srli_epi64(_a0b0h,6),11);

	a0b0 ^= _x34 ^ _x55low ^ _x21 ^ _x0;
	a0b0 &= msk;
}

inline void _reduce_94_pclmul( __m128i &a0b0 , __m128i a1b1 )
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<94>::_mask_Wbit );
	__m128i rder = _mm_load_si128( (__m128i*) gf2ext_sse<94>::_reducer_W );

	a0b0 ^= _mm_clmulepi64_si128( a1b1 , rder , 0x00 );
	__m128i a0b0h = _mm_srli_epi64( a0b0 , 30 );
	a0b0 ^= _mm_clmulepi64_si128( a0b0h , rder , 0x11 );
	a0b0 &= msk;
}

template <> inline
void _reduce<94>( __m128i &a0b0 , __m128i a1b1 )
{
//#ifdef CONFIG_FAST_PCLMULQDQ
#ifdef CONFIG_PCLMULQDQ_REDUCE
	_reduce_94_pclmul( a0b0 , a1b1 );
#else
	_reduce_94_sse( a0b0 , a1b1 );
#endif
}



////////////////////////////    gf2127   reduce   ////////////////////////////////////////


//extern const uint32_t _mask_Wbit[];
//extern const uint64_t _reducer_W[];

/// x^127 + x + 1
/// x^128 + x^2 + x


inline
void _reduce_127_shift( __m128i &a0b0 , __m128i a1b1 )
{
///__m128i _mm_set_epi32 (int i3, int i2,    int i1, int i0);
///__m128i _mm_set_epi8 (char b15, char b14,    char b13, char b12,   char b11, char b10,   char b9, char b8,   char b7, char b6,   char b5, char b4,   char b3, char b2,   char b1, char b0);

	__m128i _irr = _mm_set_epi32( 0 , 0 , 0 , 3 );
	__m128i _7byte = _mm_set_epi8(0,0,0,0, 0,0,0,0, -1,0,0,0, 0,0,0,0);

	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<127>::_mask_Wbit );

	__m128i msb = _mm_srli_si128(a0b0,15);
	__m128i _2bit = a1b1&_7byte;
	_2bit = _mm_slli_si128( _2bit , 1 );

	a0b0 ^= _mm_slli_epi64(a1b1,1)^_mm_slli_epi64(a1b1,2)
		^ _mm_srli_epi16(_2bit,7)^_mm_srli_epi16(_2bit,6)
		^(_irr&(_mm_cmplt_epi8(msb,_mm_setzero_si128()) ));
	a0b0 &= msk;
}

inline void _reduce_127_pclmul( __m128i &a0b0 , __m128i a1b1 )
{
	__m128i msk = _mm_load_si128( (__m128i*) gf2ext_sse<127>::_mask_Wbit );
	__m128i rder = _mm_load_si128( (__m128i*) gf2ext_sse<127>::_reducer_W );

	__m128i a0b0h = _mm_clmulepi64_si128( a1b1, rder , 0x11 );
	a0b0 ^= _mm_slli_si128( a0b0h , 8 );
	a0b0 ^= _mm_clmulepi64_si128( a1b1 , rder , 0x10 )
		^ (rder& _mm_cmplt_epi8(_mm_srli_si128(a0b0,15),_mm_setzero_si128()));
	a0b0 &= msk;
}

template <> inline
void _reduce<127>( __m128i &a0b0 , __m128i a1b1 )
{
//#ifdef CONFIG_FAST_PCLMULQDQ
#ifdef CONFIG_PCLMULQDQ_REDUCE
	_reduce_127_pclmul( a0b0 , a1b1 );
#else
	_reduce_127_shift( a0b0 , a1b1 );
#endif
}



////////////////////////////    gf2128   reduce   ////////////////////////////////////////

/// x^128 + x^7 + x^2 + x + 1

template <> inline
void _reduce<128>( __m128i &a0b0 , __m128i a1b1 )
{
	__m128i rder = _mm_load_si128( (__m128i*) gf2ext_u64<128>::_irrPoly );

	__m128i a0b0h = _mm_clmulepi64_si128( a1b1, rder , 0x01 );
	//a0b0 ^= _mm_clmulepi64_si128( a0b0h , rder , 0x01 );
	//a0b0 ^= _mm_shuffle_epi32( a0b0h , 0x4f );
	//a1b1 ^= _mm_shuffle_epi32( a0b0h , 0xfe );
	//a0b0 ^= _mm_shuffle_epi32( a0b0h , 0x4f );
	//a1b1 ^= (__m128i)_mm_shuffle_pd( (__m128d)a0b0h , (__m128d)_mm_setzero_si128() , 1 ); /// --> [0,a0b0_high]
	//a0b0 ^= (__m128i)_mm_shuffle_pd( (__m128d)_mm_setzero_si128() , (__m128d)a0b0h , 0 ); /// --> [a0b1_low,0]
	a1b1 ^= _mm_srli_si128( a0b0h , 8 );
	a0b0 ^= _mm_slli_si128( a0b0h , 8 );
	a0b0 ^= _mm_clmulepi64_si128( a1b1 , rder , 0x00 );
}


#endif /// defined( CONFIG_HAS_PCLMULQDQ )




/////////////////////////////////////////////////////////////////////////////


#if defined( CONFIG_PSHUFB )



template <unsigned W> struct gf2ext_sse_tbl;

template <unsigned W> inline
gf2ext_sse_tbl<W> _mul( const gf2ext_sse_tbl<W> & a , const gf2ext_sse_tbl<W> & b );


#include "blas.h"

template <unsigned W>
struct gf2ext_sse_tbl {
	__m128i v;

	typedef gf2ext_sse_tbl gf_t;
	typedef gf2ext_u64<W> gf_u64;

	gf2ext_sse_tbl() { v^=v; }
	gf2ext_sse_tbl( __m128i a ): v(a) {}
	gf2ext_sse_tbl( const gf_t & a ): v(a.v) {}
	gf2ext_sse_tbl( const gf_u64 & a ) { *((gf_u64*)this)=a; }
	explicit gf2ext_sse_tbl(bool b) { *((gf_u64*)this)=gf_u64(b); } //{ int64_t b0=b; __m128i v0=_mm_cvtsi64_si128(b0); v0=_mm_shuffle_epi32(v0,0); v =_mm_cmpgt_epi32(v0,v0^v0); }

	static const uint8_t * isomorphism;
	static const uint8_t * isomorphism_1;
	explicit gf2ext_sse_tbl( const uint8_t * x ) {
		//*((gf_u64*)this)=gf_u64(x);
		const MAT<W,W> * m = (const MAT<W,W> *) isomorphism;
		VEC<W> v1(x);
		VEC<W> v2 = m->prod( v1 );
		*((gf_u64*)this)=gf_u64( (uint8_t*)&v2 );
	}

	gf_t squ() const { return (*this)*(*this); } /// XXX:
	gf_t inv() const { return fermat_inv<gf_t,W>::exec(*this); }

	gf_t & operator *= ( const gf_t & a ) { *this = _mul( *this, a ); return *this; }
	gf_t & operator ^= ( const gf_t & a ) { v^=a.v; return *this; }
	gf_t & operator &= ( const gf_t & a ) { v&=a.v; return *this; }

	gf_t operator *( const gf_t & a ) const { return _mul(*this,a); }
	gf_t operator ^( const gf_t & a ) const { gf_t r=*this; r^=a; return r; }
	gf_t operator &( const gf_t & a ) const { gf_t r=*this; r&=a; return r; }

	gf_t & set_zero() { v^=v; return *this; }
	bool is_zero() const { return 0xffff==_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_setzero_si128(),v)); }
	//{ return ((const gf_u64 *)this)->is_zero(); }
	bool is_one() const { return ((const gf_u64 *)this)->is_one(); }

	static gf_t rand() { gf_t r = gf_u64::rand(); return r; }

	static const gf_t & zero() { return *(const gf_t* )&gf_u64::zero(); }
	static const gf_t & one() { return *(const gf_t *)&gf_u64::one(); }
	//static const gf_t & irrPoly() { return *(const gf_t *)&gf_u64::irrPoly(); }

	static gf_t assign( const uint8_t * x ) { return gf_t(x); }
	static unsigned num_byte() { return gf_u64::num_byte(); }
	void dump( uint8_t * x ) const {
		//((const gf_u64 *)this)->dump(x);
		const MAT<W,W> * m = (const MAT<W,W> *) isomorphism_1;
		const VEC<W> * v1 = (const VEC<W> *)this;
		VEC<W> v2 = m->prod( *v1 );
		v2.dump(x);
	}
	void fdump(FILE *fp) const { ((const gf_u64 *)this)->fdump(fp); }
	void fdump2(FILE *fp) const { ((const gf_u64 *)this)->fdump2(fp); }

	static const __m128i msk_low;
	static const __m128i log_tbl;
	static const __m128i exp_tbl;
	static const __m128i mul_tbl_0x2;
	static const __m128i mul_tbl_0x8;

	static const __m128i msk_5th_bit;

	static const __m128i shf_idx1;
	static const __m128i shf_idx2;
	static const __m128i shf_idx3;

	static const __m128i msk_12_bytes;
};

template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::msk_low = _mm_set_epi8(0xf,0xf,0xf,0xf, 0xf,0xf,0xf,0xf, 0xf,0xf,0xf,0xf, 0xf,0xf,0xf,0xf);

template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::log_tbl = _mm_set_epi8(12,11,13,6, 7,9,14,3, 10,5,8,2, 4,1,0,0xc0);
///{0xe0,0x0,0x1,0x4,0x2,0x8,0x5,0xa,0x3,0xe,0x9,0x7,0x6,0xd,0xb,0xc};

template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::exp_tbl = _mm_set_epi8(0x1,0x9,0xd,0xf, 0xe,0x7,0xa,0x5, 0xb,0xc,0x6,0x3, 0x8,0x4,0x2,0x1);
///{0x1,0x2,0x4,0x8, 0x3,0x6,0xc,0xb, 0x5,0xa,0x7,0xe, 0xf,0xd,0x9,0x1, 0x2,0x4,0x8,0x3,0x6,0xc,0xb,0x5,0xa,0x7,0xe,0xf,0xd,0x9,0x1,0x2};

template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::mul_tbl_0x2 = _mm_set_epi8(0xd,0xf,9,0xb, 5,7,1,3, 0xe,0xc,0xa,8, 6,4,2,0);

template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::mul_tbl_0x8 = _mm_set_epi8( 1,9,2,0xa, 7,0xf,4,0xc, 0xd,5,0xe,6, 0xb,3,8,0);


template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::msk_5th_bit = _mm_set_epi8(0x10,0x10,0x10,0x10, 0x10,0x10,0x10,0x10, 0x10,0x10,0x10,0x10, 0x10,0x10,0x10,0x10);

template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::msk_12_bytes = _mm_set_epi8(0,0,0,0, 0xff,0xff,0xff,0xff, 0xff,0xff,0xff,0xff, 0xff,0xff,0xff,0xff);

template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::shf_idx1 = _mm_set_epi8(1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1);
template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::shf_idx2 = _mm_set_epi8(2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2);
template <unsigned W>
const __m128i gf2ext_sse_tbl<W>::shf_idx3 = _mm_set_epi8(3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3);



static inline __m128i _mul_12x4( __m128i la0 , __m128i la1 , __m128i la2 , __m128i la3 , __m128i lb )
{
	typedef gf2ext_sse_tbl<96> gf_t;
	__m128i b0 = _mm_shuffle_epi8( lb , _mm_setzero_si128() );
	__m128i b1 = _mm_shuffle_epi8( lb , gf_t::shf_idx1 );
	__m128i b2 = _mm_shuffle_epi8( lb , gf_t::shf_idx2 );
	__m128i b3 = _mm_shuffle_epi8( lb , gf_t::shf_idx3 );

	b0 = _mm_add_epi8( b0 , la0 );

	b1 = _mm_add_epi8( b1 , la1 );
	b2 = _mm_add_epi8( b2 , la2 );
	b3 = _mm_add_epi8( b3 , la3 );
	b0 = _mm_add_epi8( b0 , _mm_srli_epi16( b0& gf_t::msk_5th_bit, 4) );

	b1 = _mm_add_epi8( b1 , _mm_srli_epi16( b1& gf_t::msk_5th_bit, 4) );
	b2 = _mm_add_epi8( b2 , _mm_srli_epi16( b2& gf_t::msk_5th_bit, 4) );
	b3 = _mm_add_epi8( b3 , _mm_srli_epi16( b3& gf_t::msk_5th_bit, 4) );

	__m128i r = _mm_shuffle_epi8( gf_t::exp_tbl , b0 ) ^ _mm_shuffle_epi8( gf_t::exp_tbl , b1 )
			^ _mm_shuffle_epi8( gf_t::exp_tbl , b2 ) ^_mm_shuffle_epi8( gf_t::exp_tbl , b3 );
	return r;
}

static inline void _mul_12x12( __m128i & r0, __m128i & r1 , __m128i a , __m128i b )
{
	typedef gf2ext_sse_tbl<96> gf_t;
	__m128i la = _mm_shuffle_epi8( gf_t::log_tbl , a );
	__m128i lb = _mm_shuffle_epi8( gf_t::log_tbl , b );

	__m128i la1 = _mm_alignr_epi8( la , lb , 15 );
	__m128i la2 = _mm_alignr_epi8( la , lb , 14 );
	__m128i la3 = _mm_alignr_epi8( la , lb , 13 );

	__m128i r03 = _mul_12x4( la , la1 , la2 , la3 , lb );
	lb = _mm_srli_si128( lb , 4 );
	__m128i r47 = _mul_12x4( la , la1 , la2 , la3 , lb );
	lb = _mm_srli_si128( lb , 4 );
	__m128i r8b = _mul_12x4( la , la1 , la2 , la3 , lb );

	r0 = r03 ^ _mm_slli_si128( r47 , 4 ) ^ _mm_slli_si128( r8b , 8 );
	r0 &= gf_t::msk_12_bytes;

	r1 = _mm_srli_si128( r03 , 12 ) ^ _mm_srli_si128( r47 , 8 ) ^ _mm_srli_si128( r8b , 4 );
}

static inline void _mul_12x12( __m128i & r0, __m128i & r1 , __m128i a )
{
	typedef gf2ext_sse_tbl<96> gf_t;
	__m128i la = _mm_shuffle_epi8( gf_t::log_tbl , a );
	//__m128i lb = _mm_shuffle_epi8( gf_t::log_tbl , b );
	__m128i lb = la;

	__m128i la1 = _mm_alignr_epi8( la , lb , 15 );
	__m128i la2 = _mm_alignr_epi8( la , lb , 14 );
	__m128i la3 = _mm_alignr_epi8( la , lb , 13 );

	__m128i r03 = _mul_12x4( la , la1 , la2 , la3 , lb );
	lb = _mm_srli_si128( lb , 4 );
	__m128i r47 = _mul_12x4( la , la1 , la2 , la3 , lb );
	lb = _mm_srli_si128( lb , 4 );
	__m128i r8b = _mul_12x4( la , la1 , la2 , la3 , lb );

	r0 = r03 ^ _mm_slli_si128( r47 , 4 ) ^ _mm_slli_si128( r8b , 8 );
	r0 &= gf_t::msk_12_bytes;

	r1 = _mm_srli_si128( r03 , 12 ) ^ _mm_srli_si128( r47 , 8 ) ^ _mm_srli_si128( r8b , 4 );

}



/// X^24 = a^3 X^3 + X + a
static inline void _reduce_24( __m128i & r0 , __m128i & r1 , __m128i r2 )
{
	typedef gf2ext_sse_tbl<96> gf_t;
	r0 ^= _mm_slli_si128( r2 , 1 );
	//__m128i lr2 = _mm_shuffle_epi8( gf_t::log_tbl , r2 );
	r0 ^= _mm_shuffle_epi8( gf_t::mul_tbl_0x2 , r2 );

	//__m128i lr2p3 = _mm_add_epi8( lr2 , gf_t::shf_idx3 );
	//lr2p3 = _mm_add_epi8( lr2p3 , _mm_srli_epi16( lr2p3 & gf_t::msk_5th_bit, 4) );
	//__m128i a3X3 = _mm_shuffle_epi8( gf_t::exp_tbl , lr2p3 );
	r0 ^= _mm_slli_si128( _mm_shuffle_epi8(gf_t::mul_tbl_0x8,r2) , 3);

	r1 ^= _mm_srli_si128( r0 , 12 );
	r0 &= gf_t::msk_12_bytes;
}

template <unsigned W>
gf2ext_sse_tbl<W> _mul( const gf2ext_sse_tbl<W> & a , const gf2ext_sse_tbl<W> & b )
{
	gf2ext_u64<W> r = _mul( (const gf2ext_u64<W> &) a , (const gf2ext_u64<W> &) b );
	return r;
}

template <> inline
gf2ext_sse_tbl<96> _mul( const gf2ext_sse_tbl<96> & a , const gf2ext_sse_tbl<96> & b )
{
	typedef gf2ext_sse_tbl<96> gf_t;

	__m128i a0 = a.v & gf_t::msk_low;
	__m128i a1 = _mm_srli_epi16(a.v,4) & gf_t::msk_low;
	__m128i b0 = b.v & gf_t::msk_low;
	__m128i b1 = _mm_srli_epi16(b.v,4) & gf_t::msk_low;

	__m128i r0,r1,tmp0,tmp1,rd0,rd1;
	_mul_12x12( r0, r1 , a0 , b0 );
	_mul_12x12( rd0, rd1 , a1 , b1 );
	_mul_12x12( tmp0 , tmp1 , a0^a1 , b0^b1 );

	tmp0 ^= r0 ^ rd0;
	tmp1 ^= r1 ^ rd1;
	r1 ^= tmp0;
	rd0 ^= tmp1;

	_reduce_24( r1 , rd0 , rd1 );
	_reduce_24( r0 , r1 , rd0 );

	gf_t r;
	r.v = r0|_mm_slli_epi16(r1,4);

	return r;
}

template <> inline
gf2ext_sse_tbl<96> gf2ext_sse_tbl<96>::squ() const
{
#ifdef CONFIG_PROFILE
	count_squ();
#endif
	typedef gf2ext_sse_tbl<96> gf_t;

	__m128i a0 = v & gf_t::msk_low;
	__m128i a1 = _mm_srli_epi16(v,4) & gf_t::msk_low;

	__m128i r0,r1,rd0,rd1;
	_mul_12x12( r0, r1 , a0 );
	_mul_12x12( rd0, rd1 , a1 );

	_reduce_24( r1 , rd0 , rd1 );
	_reduce_24( r0 , r1 , rd0 );

	gf_t r;
	r.v = r0|_mm_slli_epi16(r1,4);

	return r;

}

/////////////////////////////  gf 2^4^32  ===   gf 2^128  ////////////////////////////

static inline __m128i _mul_8x8( __m128i la , __m128i lb )
{
	typedef gf2ext_sse_tbl<128> gf_t;
	__m128i idx0 = la^la;
	__m128i idx4 = _mm_add_epi8(gf_t::shf_idx2,gf_t::shf_idx2);
	__m128i idx5 = _mm_add_epi8(idx4,gf_t::shf_idx1);
	__m128i idx6 = _mm_add_epi8(idx4,gf_t::shf_idx2);
	__m128i idx7 = _mm_add_epi8(idx4,gf_t::shf_idx3);

	__m128i tmp0 = _mm_add_epi8(la,_mm_shuffle_epi8(lb,idx0));
	__m128i tmp1 = _mm_add_epi8(_mm_alignr_epi8(la,la,15),_mm_shuffle_epi8(lb,gf_t::shf_idx1));
	__m128i tmp2 = _mm_add_epi8(_mm_alignr_epi8(la,la,14),_mm_shuffle_epi8(lb,gf_t::shf_idx2));
	__m128i tmp3 = _mm_add_epi8(_mm_alignr_epi8(la,la,13),_mm_shuffle_epi8(lb,gf_t::shf_idx3));
	__m128i tmp4 = _mm_add_epi8(_mm_alignr_epi8(la,la,12),_mm_shuffle_epi8(lb,idx4));
	__m128i tmp5 = _mm_add_epi8(_mm_alignr_epi8(la,la,11),_mm_shuffle_epi8(lb,idx5));
	__m128i tmp6 = _mm_add_epi8(_mm_alignr_epi8(la,la,10),_mm_shuffle_epi8(lb,idx6));
	__m128i tmp7 = _mm_add_epi8(_mm_alignr_epi8(la,la,9),_mm_shuffle_epi8(lb,idx7));

	__m128i r = _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp0 , _mm_srli_epi16( tmp0& gf_t::msk_5th_bit, 4) ) )
		^ _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp1 , _mm_srli_epi16( tmp1& gf_t::msk_5th_bit, 4) ) )
		^ _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp2 , _mm_srli_epi16( tmp2& gf_t::msk_5th_bit, 4) ) )
		^ _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp3 , _mm_srli_epi16( tmp3& gf_t::msk_5th_bit, 4) ) )
		^ _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp4 , _mm_srli_epi16( tmp4& gf_t::msk_5th_bit, 4) ) )
		^ _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp5 , _mm_srli_epi16( tmp5& gf_t::msk_5th_bit, 4) ) )
		^ _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp6 , _mm_srli_epi16( tmp6& gf_t::msk_5th_bit, 4) ) )
		^ _mm_shuffle_epi8( gf_t::exp_tbl, _mm_add_epi8( tmp7 , _mm_srli_epi16( tmp7& gf_t::msk_5th_bit, 4) ) );
	return r;
}



#ifdef _DEBUG_
extern const uint8_t * _gf16_e0;
extern const uint8_t * _gf16_l0;
static inline uint8_t gf16mul( uint8_t a , uint8_t b )
{
	if( (a==0)||(b==0) ) return 0;
	return _gf16_e0[ _gf16_l0[a]+_gf16_l0[b] ];
}
static inline __m128i _mul_8x8l_ref( __m128i la , __m128i lb )
{
	uint8_t vla[16];
	uint8_t vlb[16];
	const uint8_t * ptr = (const uint8_t *)&la;
	for(unsigned i=0;i<8;i++) vla[i] = ptr[i];
	ptr = (const uint8_t *)&lb;
	for(unsigned i=0;i<8;i++) vlb[i] = ptr[i];
	uint8_t mm[16] = {0};
	for(unsigned i=0;i<8;i++) {
		for(unsigned j=0;j<8;j++) {
			if( vla[i]>=15 || vlb[j]>=15 ) continue;
			mm[i+j] ^= _gf16_e0[ vla[i]+vlb[j] ];
		}
	}
	__m128i r0;
	uint8_t * wptr = (uint8_t *) & r0;
	for(unsigned i=0;i<16;i++) wptr[i] = mm[i];
	return r0;
}
static inline __m128i _mul_8x8_ref( __m128i a , __m128i b )
{
	uint8_t va[16];
	uint8_t vb[16];
	const uint8_t * ptr = (const uint8_t *)&a;
	for(unsigned i=0;i<8;i++) va[i] = ptr[i];
	ptr = (const uint8_t *)&b;
	for(unsigned i=0;i<8;i++) vb[i] = ptr[i];
	uint8_t mm[16] = {0};
	for(unsigned i=0;i<8;i++) {
		for(unsigned j=0;j<8;j++) {
			mm[i+j] ^= gf16mul(va[i],vb[j]);
		}
	}
	__m128i r0;
	uint8_t * wptr = (uint8_t *) & r0;
	for(unsigned i=0;i<16;i++) wptr[i] = mm[i];
	return r0;
}
static inline void _mul_16x16_ref( __m128i & r0, __m128i & r1 , __m128i a , __m128i b )
{
	uint8_t va[16];
	uint8_t vb[16];

	const uint8_t * ptr = (const uint8_t *)&a;
	for(unsigned i=0;i<16;i++) va[i] = ptr[i]&0xf;
	ptr = (const uint8_t *)&b;
	for(unsigned i=0;i<16;i++) vb[i] = ptr[i]&0xf;

	uint8_t mm[32] = {0};
	for(unsigned i=0;i<16;i++) {
		for(unsigned j=0;j<16;j++) {
			mm[i+j] ^= gf16mul(va[i],vb[j]);
		}
	}
	uint8_t * wptr = (uint8_t *) & r0;
	for(unsigned i=0;i<16;i++) wptr[i] = mm[i];
	wptr = (uint8_t *) & r1;
	for(unsigned i=0;i<16;i++) wptr[i] = mm[16+i];
}
static inline bool eq( __m128i a , __m128i b )
{
	__m128i c = a^b;
	return 0xffff==_mm_movemask_epi8(_mm_cmpeq_epi8(c,a^a));
}

static inline void dump16( const char * s , __m128i a )
{
	const uint32_t * ptr = (const uint32_t *)& a;
	printf("%s, %x %x %x %x\n", s , ptr[3], ptr[2] , ptr[1] , ptr[0] );
}
#endif

static inline void _mul_16x16( __m128i & r0, __m128i & r1 , __m128i a , __m128i b )
{
	typedef gf2ext_sse_tbl<128> gf_t;
	__m128i a0 = _mm_move_epi64(a);
	__m128i b0 = _mm_move_epi64(b);
	__m128i a1 = _mm_srli_si128(a,8);
	__m128i b1 = _mm_srli_si128(b,8);

	__m128i la0 = _mm_shuffle_epi8( gf_t::log_tbl , a0 );
	__m128i la1 = _mm_shuffle_epi8( gf_t::log_tbl , a1 );
	__m128i lb0 = _mm_shuffle_epi8( gf_t::log_tbl , b0 );
	__m128i lb1 = _mm_shuffle_epi8( gf_t::log_tbl , b1 );

	__m128i a0b0 = _mul_8x8( la0 , lb0 );
	__m128i a1b1 = _mul_8x8( la1 , lb1 );

	__m128i la0pa1 = _mm_shuffle_epi8( gf_t::log_tbl , a0^a1 );
	__m128i lb0pb1 = _mm_shuffle_epi8( gf_t::log_tbl , b0^b1 );
	__m128i a0b1pa1b0 = _mul_8x8( la0pa1 , lb0pb1 ) ^ a0b0 ^ a1b1;

	r0 = a0b0 ^ _mm_slli_si128(a0b1pa1b0, 8 );
	r1 = a1b1 ^ _mm_srli_si128(a0b1pa1b0, 8 );
}

static inline void _mul_16x16( __m128i & r0, __m128i & r1 , __m128i a )
{
	typedef gf2ext_sse_tbl<128> gf_t;
	__m128i a0 = _mm_move_epi64(a);
	__m128i a1 = _mm_srli_si128(a,8);

	__m128i la0 = _mm_shuffle_epi8( gf_t::log_tbl , a0 );
	__m128i la1 = _mm_shuffle_epi8( gf_t::log_tbl , a1 );

	__m128i a0b0 = _mul_8x8( la0 , la0 );
	__m128i a1b1 = _mul_8x8( la1 , la1 );

	__m128i la0pa1 = _mm_shuffle_epi8( gf_t::log_tbl , a0^a1 );
	__m128i a0b1pa1b0 = _mul_8x8( la0pa1 , la0pa1 ) ^ a0b0 ^ a1b1;

	r0 = a0b0 ^ _mm_slli_si128(a0b1pa1b0, 8 );
	r1 = a1b1 ^ _mm_srli_si128(a0b1pa1b0, 8 );
}

/// X^32 = X^3 + X + a
static inline void _reduce_32( __m128i & r0 , __m128i & r1 , __m128i r2 )
{
	typedef gf2ext_sse_tbl<128> gf_t;
	r0 ^= _mm_slli_si128( r2 , 1 ) ^ _mm_slli_si128( r2, 3) ^ _mm_shuffle_epi8( gf_t::mul_tbl_0x2 , r2 );
	r1 ^= _mm_srli_si128( r2 , 13 ) ^ _mm_srli_si128( r2 , 15 );
}


template <> inline
gf2ext_sse_tbl<128> _mul( const gf2ext_sse_tbl<128> & a , const gf2ext_sse_tbl<128> & b )
{
	typedef gf2ext_sse_tbl<128> gf_t;

	__m128i a0 = a.v & gf_t::msk_low;
	__m128i a1 = _mm_srli_epi16(a.v,4) & gf_t::msk_low;
	__m128i b0 = b.v & gf_t::msk_low;
	__m128i b1 = _mm_srli_epi16(b.v,4) & gf_t::msk_low;

	__m128i r0,r1,tmp0,tmp1,rd0,rd1;

	_mul_16x16( r0, r1 , a0 , b0 );
	_mul_16x16( rd0, rd1 , a1 , b1 );
	_mul_16x16( tmp0 , tmp1 , a0^a1 , b0^b1 );

	tmp0 ^= r0 ^ rd0;
	tmp1 ^= r1 ^ rd1;
	r1 ^= tmp0;
	rd0 ^= tmp1;

	_reduce_32( r1 , rd0 , rd1 );
	_reduce_32( r0 , r1 , rd0 );

	gf_t r;
	r.v = r0|_mm_slli_epi16(r1,4);

	return r;
}


template <> inline
gf2ext_sse_tbl<128> gf2ext_sse_tbl<128>::squ() const
{
	typedef gf2ext_sse_tbl<128> gf_t;

	__m128i a0 = v & gf_t::msk_low;
	__m128i a1 = _mm_srli_epi16(v,4) & gf_t::msk_low;

	__m128i r0,r1,rd0,rd1;

	_mul_16x16( r0, r1 , a0 );
	_mul_16x16( rd0, rd1 , a1 );

	_reduce_32( r1 , rd0 , rd1 );
	_reduce_32( r0 , r1 , rd0 );

	gf_t r;
	r.v = r0|_mm_slli_epi16(r1,4);

	return r;
}


#endif /// defined( CONFIG_PSHUFB )




#endif  /// _GF2EXT_SSE_HPP_s
