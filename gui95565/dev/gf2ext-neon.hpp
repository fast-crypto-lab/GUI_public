#ifndef _GF2EXT_NEON_HPP_
#define _GF2EXT_NEON_HPP_


#include "arm_neon.h"

#include <stdint.h>

#include "gf2ext.hpp"

extern const uint8_t * _gf16_log_tbl;
extern const uint8_t * _gf16_exp_tbl;
extern const uint8_t * _msk_12_byte;

extern const uint8_t * _gf16_mul_0x2;
extern const uint8_t * _gf16_mul_0x8;


#ifdef _DEBUG_
extern const uint8_t * _gf16_l0;
extern const uint8_t * _gf16_e0;
#endif

template <unsigned W> struct gf2ext_neon;

template <unsigned W>
static inline
gf2ext_neon<W> _mul( const gf2ext_neon<W> & a , const gf2ext_neon<W> & b );

#include "blas.h"

template <unsigned W>
struct gf2ext_neon {
	static const uint8_t _mask[16];
	uint8x16_t v;

	typedef gf2ext_neon gf_t;
	typedef gf2ext_u64<W> gf_u64;

	gf2ext_neon() { v = veorq_u8(v,v); }
	gf2ext_neon( uint8x16_t a ): v(a) {}
	gf2ext_neon( const gf_t & a ): v(a.v) {}
	gf2ext_neon( const gf_u64 & a ) { v= vandq_u8( vld1q_u8((const uint8_t *)&a),vld1q_u8(_mask) ); }
	explicit gf2ext_neon( bool b ) { uint64_t b0[2]={0}; b0[0]-=(uint64_t)b; b0[1]=b0[0]; v=vld1q_u8((const uint8_t *)b0); }

	static const uint8_t * isomorphism;
	static const uint8_t * isomorphism_1;
	explicit gf2ext_neon( const uint8_t * x ) {
		//v= vandq_u8( vld1q_u8(x) , vld1q_u8(_mask) );
		const MAT<W,W> * m = (const MAT<W,W> *) isomorphism;
		VEC<W> v1(x);
		VEC<W> v2 = m->prod( v2 );
		v = vandq_u8( vld1q_u8((const uint8_t*)&v2) , vld1q_u8(_mask) );
	}


	gf_t squ() const { return (*this)*(*this); } /// XXX:
	gf_t inv() const { return fermat_inv<gf_t,W>::exec(*this); }

	gf_t & operator *= ( const gf_t & a ) { *this = _mul( *this, a ); return *this; }
	gf_t & operator ^= ( const gf_t & a ) { v^=a.v; return *this; }
	gf_t & operator &= ( const gf_t & a ) { v&=a.v; return *this; }

	gf_t operator *( const gf_t & a ) const { return _mul(*this,a); }
	gf_t operator ^( const gf_t & a ) const { gf_t r=*this; r^=a; return r; }
	gf_t operator &( const gf_t & a ) const { gf_t r=*this; r&=a; return r; }

	gf_t & set_zero() { v = veorq_u8(v,v); return *this; }
	bool is_zero() const { gf_u64 vu64; vst1q_u8((uint8_t*)&vu64,v); return vu64.is_zero(); }
	bool is_one() const { gf_u64 vu64; vst1q_u8((uint8_t*)&vu64,v); return vu64.is_one(); }

	static gf_t rand() { gf_t r = gf_u64::rand(); return r; }

	static const gf_t & zero() { return *(const gf_t *)&gf_u64::zero(); }
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

};

//template <unsigned W>
//const uint8_t gf2ext_neon<W>::_mask[16];





static inline uint8x16_t vtbl2q_u8( uint8x8x2_t ltb, uint8x16_t a )
{
	uint8x8_t low,high,rl,rh;
	low  = vget_low_u8( a );
	high = vget_high_u8( a );
	rl = vtbl2_u8( ltb , low );
	rh = vtbl2_u8( ltb , high );
	return vcombine_u8( rl , rh );
}

static inline uint8x16_t vtbl4q_u8( uint8x8x4_t etb, uint8x16_t a )
{
	uint8x8_t low,high,rl,rh;
	low  = vget_low_u8( a );
	high = vget_high_u8( a );
	rl = vtbl4_u8( etb , low );
	rh = vtbl4_u8( etb , high );
	return vcombine_u8( rl , rh );
}

#ifdef _DEBUG_
static inline uint8_t gf16mul( uint8_t a , uint8_t b )
{
	if( (a==0)||(b==0) ) return 0;
	return _gf16_e0[ _gf16_l0[a]+_gf16_l0[b] ];
}
static inline uint8x16_t _mul_ref( uint8x16_t a , uint8x16_t b )
{
	uint8_t va[24];
	uint8_t vb[24];

	const uint8_t * ptr = (const uint8_t *)&a;
	for(unsigned i=0;i<12;i++) va[i] = ptr[i]&0xf;
	for(unsigned i=0;i<12;i++) va[12+i] = (ptr[i]>>4)&0xf;

	ptr = (const uint8_t *)&b;
	for(unsigned i=0;i<12;i++) vb[i] = ptr[i]&0xf;
	for(unsigned i=0;i<12;i++) vb[12+i] = (ptr[i]>>4)&0xf;

	uint8_t mm[48] = {0};
	for(unsigned i=0;i<24;i++) {
		for(unsigned j=0;j<24;j++) {
			mm[i+j] ^= gf16mul(va[i],vb[j]);
		}
	}
/// X^24 = a^3 X^3 + X + a
	for(unsigned i=46;i>=24;i--){
		mm[i-24] ^= gf16mul(mm[i],2);
		mm[i-24+1] ^= mm[i];
		mm[i-24+3] ^= gf16mul(mm[i],8);
	}

	uint8_t vc[32] = {0};
	for(unsigned i=0;i<12;i++) vc[i]=mm[i];
	for(unsigned i=0;i<12;i++) vc[i] ^= (mm[12+i]<<4);
	return vld1q_u8( vc );
}
static inline uint8x16_t _mul_32x32_ref( uint8x16_t a , uint8x16_t b )
{
	uint8_t va[32];
	uint8_t vb[32];
	const uint8_t * ptr = (const uint8_t *)&a;
	for(unsigned i=0;i<16;i++) va[i] = ptr[i]&0xf;
	for(unsigned i=0;i<16;i++) va[16+i] = (ptr[i]>>4)&0xf;
	ptr = (const uint8_t *)&b;
	for(unsigned i=0;i<16;i++) vb[i] = ptr[i]&0xf;
	for(unsigned i=0;i<16;i++) vb[16+i] = (ptr[i]>>4)&0xf;
	uint8_t mm[64] = {0};
	for(unsigned i=0;i<32;i++) {
		for(unsigned j=0;j<32;j++) {
			mm[i+j] ^= gf16mul(va[i],vb[j]);
		}
	}
/// X^32 = X^3 + X + a
	for(unsigned i=63;i>=32;i--){
		mm[i-32] ^= gf16mul(mm[i],2);
		mm[i-32+1] ^= mm[i];
		mm[i-32+3] ^= mm[i];
	}
	uint8_t vc[32] = {0};
	for(unsigned i=0;i<16;i++) vc[i]=mm[i];
	for(unsigned i=0;i<16;i++) vc[i] ^= (mm[16+i]<<4);
	return vld1q_u8( vc );
}
static inline void _mul_16x16_ref( uint8x16_t & r0, uint8x16_t & r1 , uint8x16_t a , uint8x16_t b )
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
	r0 = vld1q_u8(mm);
	r1 = vld1q_u8(mm+16);
}
#endif /// def _DEBUG_




static inline uint8x16_t _mul_12x4( uint8x16_t la0 , uint8x16_t la1 , uint8x16_t la2 , uint8x16_t la3
	, uint8x16_t lb0 , uint8x16_t lb1 , uint8x16_t lb2 , uint8x16_t lb3 , uint8x8x4_t exp_tb )
{
//	uint8x8x4_t exp_tb = vld4_u8( _gf16_exp_tbl );

	lb0 = vtbl4q_u8( exp_tb , vaddq_u8( lb0 , la0 ) );
	lb1 = vtbl4q_u8( exp_tb , vaddq_u8( lb1 , la1 ) );
	lb2 = vtbl4q_u8( exp_tb , vaddq_u8( lb2 , la2 ) );
	lb3 = vtbl4q_u8( exp_tb , vaddq_u8( lb3 , la3 ) );

	return veorq_u8( veorq_u8(lb0,lb1) , veorq_u8(lb2,lb3) );
}

static inline void _mul_12x12( uint8x16_t & r0, uint8x16_t & r1 , uint8x16_t a , uint8x16_t b )
{
	uint8x8x2_t log_tb = vld2_u8( _gf16_log_tbl );
	uint8x8x4_t exp_tb = vld4_u8( _gf16_exp_tbl );
	uint8x16_t msk_12b = vld1q_u8( _msk_12_byte );

	uint8x16_t la = vtbl2q_u8( log_tb , a );

	uint8x8_t low = vget_low_u8( b );
	uint8x8_t high = vget_high_u8( b );
	uint8x8_t rl = vtbl2_u8( log_tb , low );
	uint8x8_t rh = vtbl2_u8( log_tb , high );
	//uint8x16_t lb = vcombine_u8( rl , rh );

	uint8x16_t la1 = vextq_u8(la,la,15);
	uint8x16_t la2 = vextq_u8(la,la,14);
	uint8x16_t la3 = vextq_u8(la,la,13);

	uint8x16_t lb0 = vdupq_lane_u8(rl,0);
	uint8x16_t lb1 = vdupq_lane_u8(rl,1);
	uint8x16_t lb2 = vdupq_lane_u8(rl,2);
	uint8x16_t lb3 = vdupq_lane_u8(rl,3);
	uint8x16_t r03 = _mul_12x4(la,la1,la2,la3,lb0,lb1,lb2,lb3,exp_tb);

	lb0 = vdupq_lane_u8(rl,4);
	lb1 = vdupq_lane_u8(rl,5);
	lb2 = vdupq_lane_u8(rl,6);
	lb3 = vdupq_lane_u8(rl,7);
	uint8x16_t r47 = _mul_12x4(la,la1,la2,la3,lb0,lb1,lb2,lb3,exp_tb);

	lb0 = vdupq_lane_u8(rh,0);
	lb1 = vdupq_lane_u8(rh,1);
	lb2 = vdupq_lane_u8(rh,2);
	lb3 = vdupq_lane_u8(rh,3);
	uint8x16_t r8b = _mul_12x4(la,la1,la2,la3,lb0,lb1,lb2,lb3,exp_tb);

	uint8x16_t zero = veorq_u8( a , a );
	r0 = veorq_u8( r03 , vextq_u8(zero,r47,12) );
	r0 = vandq_u8( msk_12b , veorq_u8( r0 , vextq_u8(zero,r8b,8) ) );

	r1 = veorq_u8( vextq_u8(r03,zero,12) , vextq_u8(r47,zero,8) );
	r1 = veorq_u8( r1 , vextq_u8(r8b,zero,4) );
}


static inline void _mul_12x12( uint8x16_t & r0, uint8x16_t & r1 , uint8x16_t a )
{
	uint8x8x2_t log_tb = vld2_u8( _gf16_log_tbl );
	uint8x8x4_t exp_tb = vld4_u8( _gf16_exp_tbl );
	uint8x16_t msk_12b = vld1q_u8( _msk_12_byte );

	//uint8x16_t la = vtbl2q_u8( log_tb , a );

	uint8x8_t low = vget_low_u8( a );
	uint8x8_t high = vget_high_u8( a );
	uint8x8_t rl = vtbl2_u8( log_tb , low );
	uint8x8_t rh = vtbl2_u8( log_tb , high );
	uint8x16_t la = vcombine_u8( rl , rh );

	uint8x16_t la1 = vextq_u8(la,la,15);
	uint8x16_t la2 = vextq_u8(la,la,14);
	uint8x16_t la3 = vextq_u8(la,la,13);

	uint8x16_t lb0 = vdupq_lane_u8(rl,0);
	uint8x16_t lb1 = vdupq_lane_u8(rl,1);
	uint8x16_t lb2 = vdupq_lane_u8(rl,2);
	uint8x16_t lb3 = vdupq_lane_u8(rl,3);
	uint8x16_t r03 = _mul_12x4(la,la1,la2,la3,lb0,lb1,lb2,lb3,exp_tb);

	lb0 = vdupq_lane_u8(rl,4);
	lb1 = vdupq_lane_u8(rl,5);
	lb2 = vdupq_lane_u8(rl,6);
	lb3 = vdupq_lane_u8(rl,7);
	uint8x16_t r47 = _mul_12x4(la,la1,la2,la3,lb0,lb1,lb2,lb3,exp_tb);

	lb0 = vdupq_lane_u8(rh,0);
	lb1 = vdupq_lane_u8(rh,1);
	lb2 = vdupq_lane_u8(rh,2);
	lb3 = vdupq_lane_u8(rh,3);
	uint8x16_t r8b = _mul_12x4(la,la1,la2,la3,lb0,lb1,lb2,lb3,exp_tb);

	uint8x16_t zero = veorq_u8( a , a );
	r0 = veorq_u8( r03 , vextq_u8(zero,r47,12) );
	r0 = vandq_u8( msk_12b , veorq_u8( r0 , vextq_u8(zero,r8b,8) ) );

	r1 = veorq_u8( vextq_u8(r03,zero,12) , vextq_u8(r47,zero,8) );
	r1 = veorq_u8( r1 , vextq_u8(r8b,zero,4) );
}



/// X^24 = a^3 X^3 + X + a
static inline void _reduce_24( uint8x16_t & r0 , uint8x16_t & r1 , uint8x16_t r2 )
{
	uint8x16_t zero = veorq_u8(r2,r2);
	uint8x16_t msk_12b = vld1q_u8( _msk_12_byte );
	uint8x8x2_t mul_0x2 = vld2_u8( _gf16_mul_0x2 );
	uint8x8x2_t mul_0x8 = vld2_u8( _gf16_mul_0x8 );

	uint8x16_t r2p1 = vtbl2q_u8( mul_0x2 , r2);
	uint8x16_t r2p3 = vtbl2q_u8( mul_0x8 , r2);

	r0 = veorq_u8( veorq_u8( r0,vextq_u8(zero,r2p3,13) ) , veorq_u8( vextq_u8(zero,r2,15),r2p1 ) );
	r1 = veorq_u8( r1 , vextq_u8(r0,zero,12) );
	r0 = vandq_u8( r0 , msk_12b );
}

template <unsigned W>
static inline
gf2ext_neon<W> _mul( const gf2ext_neon<W> & a , const gf2ext_neon<W> & b )
{
	gf2ext_u64<W> r = _mul( (const gf2ext_u64<W> &) a , (const gf2ext_u64<W> &) b );
	return r;
}



template <>
inline
gf2ext_neon<96> _mul( const gf2ext_neon<96> & a , const gf2ext_neon<96> & b )
{
	uint8x16_t _bit_0xf = vdupq_n_u8( 0xf );
	uint8x16_t a0 = vandq_u8( a.v , _bit_0xf );
	uint8x16_t a1 = vshrq_n_u8(a.v,4);
	uint8x16_t b0 = vandq_u8( b.v , _bit_0xf );
	uint8x16_t b1 = vshrq_n_u8(b.v,4);

	uint8x16_t r0,r1,tmp0,tmp1,rd0,rd1;
	_mul_12x12( r0, r1 , a0 , b0 );
	_mul_12x12( rd0, rd1 , a1 , b1 );
	_mul_12x12( tmp0 , tmp1 , veorq_u8(a0,a1) , veorq_u8(b0,b1) );
	tmp0 = veorq_u8( tmp0 , veorq_u8( r0 , rd0 ) );
	tmp1 = veorq_u8( tmp1 , veorq_u8( r1 , rd1 ) );
	r1 = veorq_u8( r1,tmp0 );
	rd0 = veorq_u8( rd0,tmp1 );

	_reduce_24( r1 , rd0 , rd1 );
	_reduce_24( r0 , r1 , rd0 );

	uint8x16_t r = vorrq_u8(r0, vshlq_n_u8(r1,4) );
	return gf2ext_neon<96>(r);
}



template <> inline
gf2ext_neon<96> gf2ext_neon<96>::squ() const
{
	uint8x16_t _bit_0xf = vdupq_n_u8( 0xf );
	uint8x16_t a0 = vandq_u8( v , _bit_0xf );
	uint8x16_t a1 = vshrq_n_u8(v,4);

	uint8x16_t r0,r1,rd0,rd1;
	_mul_12x12( r0, r1 , a0 );
	_mul_12x12( rd0, rd1 , a1 );

	_reduce_24( r1 , rd0 , rd1 );
	_reduce_24( r0 , r1 , rd0 );

	uint8x16_t r = veorq_u8(r0, vshlq_n_u8(r1,4) );
	return gf2ext_neon<96>(r);

}


////////////////////////////////  gf 2^128 == gf 2^4^32  ///////////////////////////////////////

static inline uint8x16_t _mul_8x8( uint8x8_t a , uint8x8_t b )
{
	uint8x8x2_t log_tb = vld2_u8( _gf16_log_tbl );
	uint8x8x4_t exp_tb = vld4_u8( _gf16_exp_tbl );

	uint8x8_t lb = vtbl2_u8( log_tb , b );
	uint8x8_t la = vtbl2_u8( log_tb , a );

	uint8x8_t r0 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,0)) );
	uint8x8_t r1 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,1)) );
	uint8x8_t r2 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,2)) );
	uint8x8_t r3 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,3)) );
	uint8x8_t r4 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,4)) );
	uint8x8_t r5 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,5)) );
	uint8x8_t r6 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,6)) );
	uint8x8_t r7 = vtbl4_u8( exp_tb , vadd_u8(la,vdup_lane_u8(lb,7)) );
	uint8x8_t zero = veor_u8( r0 , r0 );
	uint8x16_t zero2 = vcombine_u8( zero , zero );

	uint8x16_t r1_ = vextq_u8( zero2 , vcombine_u8(r1,zero) , 15 );
	uint8x16_t r2_ = vextq_u8( zero2 , vcombine_u8(r2,zero) , 14 );
	uint8x16_t r3_ = vextq_u8( zero2 , vcombine_u8(r3,zero) , 13 );
	uint8x16_t r4_ = vextq_u8( zero2 , vcombine_u8(r4,zero) , 12 );
	uint8x16_t r5_ = vextq_u8( zero2 , vcombine_u8(r5,zero) , 11 );
	uint8x16_t r6_ = vextq_u8( zero2 , vcombine_u8(r6,zero) , 10 );
	uint8x16_t r7_ = vextq_u8( zero2 , vcombine_u8(r7,zero) , 9 );

	return veorq_u8( veorq_u8( veorq_u8( vcombine_u8(r0,zero),r1_) , veorq_u8(r2_,r3_) )
			, veorq_u8( veorq_u8(r4_,r5_) , veorq_u8(r6_,r7_) ) );
}

static inline void _mul_16x16( uint8x16_t & r0, uint8x16_t & r1 , uint8x16_t a , uint8x16_t b )
{

	uint8x8_t a0 = vget_low_u8( a );
	uint8x8_t a1 = vget_high_u8( a );
	uint8x8_t b0 = vget_low_u8( b );
	uint8x8_t b1 = vget_high_u8( b );
	uint8x8_t zero = veor_u8( a0,a0 );

	uint8x16_t a0b0 = _mul_8x8( a0 , b0 );
	uint8x16_t a1b1 = _mul_8x8( a1 , b1 );

	uint8x16_t a0b1pa1b0 = _mul_8x8( veor_u8(a0,a1), veor_u8(b0,b1) ) ^ a0b0 ^ a1b1 ;
	r0 = veorq_u8( a0b0 , vcombine_u8( zero , vget_low_u8(a0b1pa1b0) ) );
	r1 = veorq_u8( a1b1 , vcombine_u8( vget_high_u8(a0b1pa1b0) , zero ) );
}

static inline void _mul_16x16( uint8x16_t & r0, uint8x16_t & r1 , uint8x16_t a )
{

	uint8x8_t a0 = vget_low_u8( a );
	uint8x8_t a1 = vget_high_u8( a );
	uint8x8_t zero = veor_u8( a0,a0 );

	uint8x16_t a0a0 = _mul_8x8( a0 , a0 );
	uint8x16_t a1a1 = _mul_8x8( a1 , a1 );

	uint8x16_t a0b1pa1b0 = _mul_8x8( veor_u8(a0,a1), veor_u8(a0,a1) ) ^ a0a0 ^ a1a1 ;
	r0 = veorq_u8( a0a0 , vcombine_u8( zero , vget_low_u8(a0b1pa1b0) ) );
	r1 = veorq_u8( a1a1 , vcombine_u8( vget_high_u8(a0b1pa1b0) , zero ) );
}


/// X^32 = X^3 + X + a
static inline void _reduce_32( uint8x16_t & r0 , uint8x16_t & r1 , uint8x16_t r2 )
{
	uint8x16_t zero = veorq_u8(r2,r2);
	uint8x8x2_t mul_0x2 = vld2_u8( _gf16_mul_0x2 );

	r0 = veorq_u8( veorq_u8( vtbl2q_u8(mul_0x2,r2) , vextq_u8(zero,r2,15) ) , veorq_u8( vextq_u8(zero,r2,13) , r0 ) );
	r1 = veorq_u8( veorq_u8( r1 , vextq_u8(r2,zero,15)) , vextq_u8(r2,zero,13) );
}

template <>
inline
gf2ext_neon<128> _mul( const gf2ext_neon<128> & a , const gf2ext_neon<128> & b )
{

	uint8x16_t _bit_0xf = vdupq_n_u8( 0xf );
	uint8x16_t a0 = vandq_u8( a.v , _bit_0xf );
	uint8x16_t a1 = vshrq_n_u8(a.v,4);
	uint8x16_t b0 = vandq_u8( b.v , _bit_0xf );
	uint8x16_t b1 = vshrq_n_u8(b.v,4);

	uint8x16_t r0,r1,tmp0,tmp1,rd0,rd1;

	_mul_16x16( r0, r1 , a0 , b0 );
	_mul_16x16( rd0, rd1 , a1 , b1 );
	_mul_16x16( tmp0 , tmp1 , veorq_u8(a0,a1) , veorq_u8(b0,b1) );

	tmp0 = veorq_u8( tmp0 , veorq_u8( r0 , rd0 ) );
	tmp1 = veorq_u8( tmp1 , veorq_u8( r1 , rd1 ) );
	r1 = veorq_u8( r1,tmp0 );
	rd0 = veorq_u8( rd0,tmp1 );

	_reduce_32( r1 , rd0 , rd1 );
	_reduce_32( r0 , r1 , rd0 );

	uint8x16_t r = vorrq_u8(r0, vshlq_n_u8(r1,4) );

	return gf2ext_neon<128>(r);
}


template <> inline
gf2ext_neon<128> gf2ext_neon<128>::squ() const
{
	uint8x16_t _bit_0xf = vdupq_n_u8( 0xf );
	uint8x16_t a0 = vandq_u8( v , _bit_0xf );
	uint8x16_t a1 = vshrq_n_u8(v,4);

	uint8x16_t r0,r1,rd0,rd1;
	_mul_16x16( r0, r1 , a0 );
	_mul_16x16( rd0, rd1 , a1 );

	_reduce_32( r1 , rd0 , rd1 );
	_reduce_32( r0 , r1 , rd0 );

	uint8x16_t r = veorq_u8(r0, vshlq_n_u8(r1,4) );
	return gf2ext_neon<128>(r);

}




#endif  /// _GF2EXT_SSE_HPP_s
