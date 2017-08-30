#ifndef _GF2EXT_HPP_
#define _GF2EXT_HPP_

#include <stdint.h>
#include <stdio.h>

#include "crand.h"


template <unsigned W> struct gf2ext_u64;

template <unsigned W> inline
gf2ext_u64<W> _mul( const gf2ext_u64<W> & a, const gf2ext_u64<W> & b );


template <typename gf_t,unsigned W>
struct fermat_inv {
	static const gf_t exec( const gf_t & x );
};


template <unsigned W> inline
gf2ext_u64<W> _inv( const gf2ext_u64<W> & a );



template <unsigned W>
struct gf2ext_u64 {
	uint64_t v[2];

	typedef gf2ext_u64 gfext;
	typedef gf2ext_u64 gf_t;

	static const unsigned _num_byte = (W+7)/8;

	static const uint64_t _zero[2];
	static const uint64_t _one[2];
	static const uint64_t _msb_u64 __attribute__((aligned(16))) = 0x8000000000000000ULL;

	static const uint64_t _mask_high_limb __attribute__((aligned(16))) = 0xffffffffffffffffULL>>(128-W);
	static const unsigned _nbit_high_limb = (W-64);
	static const uint64_t _irrPoly[2];

	gf2ext_u64() { v[0]^=v[0]; v[1]^=v[1]; }

	explicit gf2ext_u64(bool b) { uint64_t b1=b; uint64_t b2=0; b2-=b1; v[0]=b2; v[1]=b2; v[1]&=_mask_high_limb; }

	gf2ext_u64(const gf2ext_u64 & a ) { v[0]=a.v[0]; v[1]=a.v[1]; }

	explicit gf2ext_u64( const uint8_t * x ) {
		v[0]=((const uint64_t *)x)[0];
		for(unsigned i=8;i<_num_byte;i++) ((uint8_t*)v)[i]=x[i];
		v[1] &= _mask_high_limb;
	}

	//gf_t inv() const { return _inv( *this );} /// XXX:
	gf_t inv() const { return fermat_inv<gf_t,W>::exec(*this); }
	gf_t squ() const { return _mul( *this , *this ); } /// XXX:

	gf_t & operator *= ( const gf_t & a ) { *this = _mul( *this, a ); return *this; }
	gf_t & operator ^= ( const gf_t & a ) { v[0]^=a.v[0]; v[1]^=a.v[1]; return *this; }
	gf_t & operator &= ( const gf_t & a ) { v[0]&=a.v[0]; v[1]&=a.v[1]; return *this; }

	gf_t operator *( const gf_t & a ) const { return _mul(*this,a); }
	gf_t operator ^( const gf_t & a ) const { gf_t r=*this; r^=a; return r; }
	gf_t operator &( const gf_t & a ) const { gf_t r=*this; r&=a; return r; }

	gf_t & set_zero() { v[0]^=v[0]; v[1]^=v[1]; return *this; }
	bool is_zero() const { return 0==(v[0]|v[1]); }
	bool is_one() const { return 0==((v[0]-1)|v[1]); }

	static gf_t rand() {
		gf_t r;
		RAND_bytes( (unsigned char*)r.v, 16 );
		r.v[1] &= _mask_high_limb;
		return r;
	}
	static const gf_t& zero() { return (const gf_t&)_zero; }
	static const gf_t& one() { return (const gf_t&)_one; }
	static const gf_t& irrPoly() { return (const gf_t&)_irrPoly; }

	static gf_t assign( const uint8_t * x ) { return gf_t(x); }

	static unsigned num_byte() { return _num_byte; }
	void dump( uint8_t * x ) const { for(unsigned i=0;i<_num_byte;i++) x[i]=((uint8_t*)v)[i]; }

	void fdump(FILE *fp) const {
		uint16_t *v16 = (uint16_t *)v;
		if(v[1]) fprintf(fp,".%04x.",v16[4]);
		if(v[0]) fprintf(fp,".%04x",v16[0]);
		if( is_zero() ) fprintf(fp,"0");
	}

	void fdump2(FILE *fp) const {
		uint16_t *v16 = (uint16_t *)v;
		for(int i=7;i>=0;i--) fprintf(fp,"[%4x]",v16[i]);
//		if(v[1]) fprintf(fp,".%04x.",v16[4]);
//		if(v[0]) fprintf(fp,".%04x",v16[0]);
//		if( is_zero() ) fprintf(fp,"0");
	}

//// private

	gf2ext_u64 & mul_2() { /// XXX: W = 128
		v[1] <<= 1;
		//v[1] |= (v[0]&_msb_u64)? 1:0;
		v[1] |= v[0]>>63;
		v[0] <<= 1;
		//(*this) ^= ((v[1] & _irrPoly[1])? (const gf2ext_u64 &)_irrPoly : (const gf2ext_u64 &)_zero);
		uint64_t msk = 0;
		msk -= (v[1]>>_nbit_high_limb);
		v[1] ^= _irrPoly[1]&msk;
		v[0] ^= _irrPoly[0]&msk;
		return *this;
	}

	int left_most_bit() const { /// non-time constant
		if(v[1]) return 64+63-__builtin_clzll(v[1]);
		if(v[0]) return 63-__builtin_clzll(v[0]);
		return -1;
	}

	gf2ext_u64 & bit_shl( unsigned i) { /// non-time constant
		if( i >= 64 ) { v[1] = v[0]; v[0] = 0; i -= 64; }
		v[1] = (v[1]<<i)| (v[0]>>(64-i));
		v[0] <<= i;
		return *this;
	}
};



template <unsigned W>
const uint64_t gf2ext_u64<W>::_zero[2] __attribute__((aligned(32))) = {0};

template <unsigned W>
const uint64_t gf2ext_u64<W>::_one[2] __attribute__((aligned(32))) = {1ULL,0ULL};




#include "run_config.h"

#ifdef CONFIG_PROFILE
#include "stat_profile.h"
#endif

template <> inline
gf2ext_u64<128> & gf2ext_u64<128>::mul_2()
{
	//uint64_t irr = (v[1]&0x8000000000000000ULL)? _irrPoly[0]:0;
	uint64_t irr = (0ULL - (v[1]>>63)) & _irrPoly[0];
	v[1] = (v[1]<<1) | (v[0]>>63);
	v[0] = (v[0]<<1) ^ irr;
	return *this;
}

template <unsigned W> inline
gf2ext_u64<W> _mul( const gf2ext_u64<W> & a, const gf2ext_u64<W> & b )
{
#ifdef CONFIG_PROFILE
count_mul();
#endif
	gf2ext_u64<W> bpower = b;
	gf2ext_u64<W> accu;

	uint64_t a0 = a.v[0];
	uint64_t msk = 0;
	//for( int i=0;i<64;i++) {
	//	if( a0&0x01 ) accu ^= bpower;
	for( unsigned i=64;i>0;i--){
		msk ^= msk;
		msk -= (a0&0x1);
		accu.v[1] ^= bpower.v[1]&msk;
		accu.v[0] ^= bpower.v[0]&msk;
		a0 >>= 1;
		bpower.mul_2();
	}
	a0 = a.v[1];
	//for( unsigned i=0;i< gf2ext_u64<W>::_nbit_high_limb;i++) {
	//	if( a0&0x01 ) accu ^= bpower;
	for( unsigned i=gf2ext_u64<W>::_nbit_high_limb-1; i>0 ;i--) {
		msk ^= msk;
		msk -= (a0&0x1);
		accu.v[1] ^= bpower.v[1]&msk;
		accu.v[0] ^= bpower.v[0]&msk;
		a0 >>= 1;
		bpower.mul_2();
	}
	msk ^= msk;
	msk -= (a0&0x1);
	accu.v[1] ^= bpower.v[1]&msk;
	accu.v[0] ^= bpower.v[0]&msk;
	return accu;
}


#ifdef _MUL_U64_KARATSUBA_
inline uint64_t _mul_32x32( uint64_t a , uint64_t b)
{
	uint64_t r=0;
	a &= 0xffffffffULL;
	for(unsigned i=32;i!=0;i--) {
		r ^= (b&1)? a: 0;
		b >>= 1;
		a <<= 1;
	}
	return r;
}

inline void _mul_64x64( uint64_t & rh, uint64_t & rl, uint64_t a , uint64_t b)
{
	uint64_t ah = a>>32;
	uint64_t al = a&0xffffffffULL;
	uint64_t bh = b>>32;
	uint64_t bl = b&0xffffffffULL;
	rh = _mul_32x32( ah , bh );
	rl = _mul_32x32( al , bl );
	uint64_t mid = _mul_32x32( ah^al , bh^bl ) ^ rh ^ rl;
	rl ^= (mid<<32);
	rh ^= (mid>>32);
}

/// X^128 + X^7 + X^2 + X + 1
template <> inline
gf2ext_u64<128> _mul( const gf2ext_u64<128> & a, const gf2ext_u64<128> & b )
{
	gf2ext_u64<128> r;
	_mul_64x64( r.v[1] , r.v[0] , a.v[0] , b.v[0] );
	uint64_t rh,rl,midh,midl;
	_mul_64x64(rh,rl, a.v[1] , b.v[1] );
	_mul_64x64( midh, midl , a.v[0]^a.v[1] , b.v[0]^b.v[1] );
	midh ^= r.v[1]^rh;
	midl ^= r.v[0]^rl;
	rl ^= midh;
	r.v[1] ^= midl;

	uint64_t acc = gf2ext_u64<128>::_irrPoly[0];
	uint64_t idx = 0x100000000000000ULL;
	uint64_t rd_h8 = 0;
	for(unsigned i=8;i!=0;i--) {
		rd_h8 ^= (idx&rh)? acc: 0;
		acc <<= 1;
		idx <<= 1;
	}
	r.v[1] ^= (rd_h8<<56);
	rl ^= (rd_h8>>8);

	acc = gf2ext_u64<128>::_irrPoly[0];
	idx = 0x100000000000000ULL;
	rd_h8 = 0;
	for(unsigned i=8;i!=0;i--) {
		rd_h8 ^= (idx&rl)? acc: 0;
		acc <<= 1;
		idx <<= 1;
	}
	r.v[0] ^= (rd_h8<<56);
	r.v[1] ^= (rd_h8>>8);

	acc = gf2ext_u64<128>::_irrPoly[0];
	idx = 1;
	uint64_t rd = 0;
	for(unsigned i=56;i!=0;i--) {
		rd ^= (idx&rh)? acc: 0;
		acc <<= 1;
		idx <<= 1;
	}
	r.v[1] ^= rd;

	acc = gf2ext_u64<128>::_irrPoly[0];
	idx = 1;
	rd = 0;
	for(unsigned i=56;i!=0;i--) {
		rd ^= (idx&rl)? acc: 0;
		acc <<= 1;
		idx <<= 1;
	}
	r.v[0] ^= rd;

	return r;
}
#endif

///
/// ExtGCD
///
template <unsigned W> inline
gf2ext_u64<W> _inv( const gf2ext_u64<W> & __a )
{
	typedef gf2ext_u64<W> gfext;
#ifdef CONFIG_PROFILE
count_inv();
#endif
	gfext buf1[2];
	gfext buf2[2];
	gfext * ptr1 = buf1;
	gfext * ptr2 = buf2;
	gfext * tmp;

	buf1[0] = gfext::one();
	buf1[1] = __a;
	int lmb1 = __a.left_most_bit();
	buf2[0] = gfext::zero();
	buf2[1] = gfext::irrPoly();
	int lmb2 = W;
	gfext a;
	gfext b;
	while( lmb1 > 0 ) {
		a = ptr1[0];
		b = ptr1[1];
		unsigned diff = lmb2-lmb1;
		if(diff) a.bit_shl(diff);
		if(diff) b.bit_shl(diff);
		ptr2[0] ^= a;
		ptr2[1] ^= b;
		lmb2 = ptr2[1].left_most_bit();
		if( lmb2 < lmb1 ) {
			int t = lmb1; lmb1=lmb2; lmb2=t;
			tmp = ptr1; ptr1=ptr2; ptr2=tmp;
		}
	}
	return ptr1[0];
}




template <typename gf_t>
struct fermat_inv<gf_t, 94> {
static const gf_t exec( const gf_t & x ) {
	/// Extend X to X^( 2^94 - 2 )
	gf_t x2 = x.squ();
	gf_t x3 = _mul( x2 , x );

	gf_t tmp = x3;
	for(int i=2;i>0;i--) tmp = tmp.squ();
	gf_t x_2_4_1 = _mul( tmp,x3); /// F

	tmp = x_2_4_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_8_1 = _mul( tmp,x_2_4_1); /// FF

	tmp = x_2_8_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_16_1 = _mul( tmp,x_2_8_1);  /// FFFF

	tmp = x_2_16_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_32_1 = _mul( tmp,x_2_16_1);  /// FFFF,FFFF

	tmp = x_2_32_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_64_1 = _mul( tmp,x_2_32_1);  /// FFFF,FFFF,FFFF,FFFF

	tmp = x_2_64_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_80_1 = _mul( tmp , x_2_16_1); /// FFFF,FFFF,FFFF,FFFF,FFFF

	tmp = x_2_80_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_88_1 = _mul( tmp , x_2_8_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FF

	tmp = x_2_88_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_92_1 = _mul( tmp , x_2_4_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FFF

	tmp = x_2_92_1;
	for(int i=2;i>0;i--) tmp = tmp.squ();
	return _mul( tmp , x2 );
}
};


template <typename gf_t>
struct fermat_inv<gf_t, 95> {
static const gf_t exec( const gf_t & x ) {
	/// Extend X to X^( 2^95 - 2 )
	gf_t x2 = x.squ();
	gf_t x3 = _mul( x2 , x );

	gf_t tmp = x3;
	tmp = tmp.squ();
	gf_t x_2_3_2 = tmp; ///
	tmp = tmp.squ();
	gf_t x_2_4_1 = _mul( tmp,x3); /// F

	tmp = x_2_4_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_8_1 = _mul( tmp,x_2_4_1); /// FF

	tmp = x_2_8_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_16_1 = _mul( tmp,x_2_8_1);  /// FFFF

	tmp = x_2_16_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_32_1 = _mul( tmp,x_2_16_1);  /// FFFF,FFFF

	tmp = x_2_32_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_64_1 = _mul( tmp,x_2_32_1);  /// FFFF,FFFF,FFFF,FFFF

	tmp = x_2_64_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_80_1 = _mul( tmp , x_2_16_1); /// FFFF,FFFF,FFFF,FFFF,FFFF

	tmp = x_2_80_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_88_1 = _mul( tmp , x_2_8_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FF

	tmp = x_2_88_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_92_1 = _mul( tmp , x_2_4_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FFF

	tmp = x_2_92_1;
	for(int i=3;i>0;i--) tmp = tmp.squ();
	return _mul( tmp , x_2_3_2 );
}
};



template <typename gf_t>
struct fermat_inv<gf_t, 96> {
static const gf_t exec( const gf_t & x ) {
	/// Extend X to X^( 2^96 - 2 )
	gf_t x2 = x.squ();
	gf_t x3 = _mul( x2 , x );

	gf_t tmp = x3;
	for(int i=2;i>0;i--) tmp = tmp.squ();
	gf_t x_2_4_1 = _mul( tmp,x3); /// F
	gf_t x_2_4_2 = _mul( tmp,x2); /// E

	tmp = x_2_4_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_8_1 = _mul( tmp,x_2_4_1); /// FF

	tmp = x_2_8_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_16_1 = _mul( tmp,x_2_8_1);  /// FFFF

	tmp = x_2_16_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_32_1 = _mul( tmp,x_2_16_1);  /// FFFF,FFFF

	tmp = x_2_32_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_64_1 = _mul( tmp,x_2_32_1);  /// FFFF,FFFF,FFFF,FFFF

	tmp = x_2_64_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_80_1 = _mul( tmp , x_2_16_1); /// FFFF,FFFF,FFFF,FFFF,FFFF

	tmp = x_2_80_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_88_1 = _mul( tmp , x_2_8_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FF

	tmp = x_2_88_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_92_1 = _mul( tmp , x_2_4_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FFF

	tmp = x_2_92_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	return _mul( tmp , x_2_4_2 );
}
};


template <typename gf_t>
struct fermat_inv<gf_t, 103> {
static const gf_t exec( const gf_t & x ) {
	/// Extend X to X^( 2^103 - 2 )
	gf_t x2 = x.squ();
	gf_t x3 = _mul( x2 , x );

	gf_t tmp = x3;
	tmp = tmp.squ();
	gf_t x_2_3_2 = tmp; ///
	tmp = tmp.squ();
	gf_t x_2_4_1 = _mul( tmp,x3); /// F

	tmp = x_2_4_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_8_1 = _mul( tmp,x_2_4_1); /// FF

	tmp = x_2_8_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_16_1 = _mul( tmp,x_2_8_1);  /// FFFF

	tmp = x_2_16_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_32_1 = _mul( tmp,x_2_16_1);  /// FFFF,FFFF

	tmp = x_2_32_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_64_1 = _mul( tmp,x_2_32_1);  /// FFFF,FFFF,FFFF,FFFF

	tmp = x_2_64_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_96_1 = _mul( tmp , x_2_32_1); /// FFFFx6

	tmp = x_2_96_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_100_1 = _mul( tmp , x_2_4_1); /// FFFFx6,F

	tmp = x_2_100_1;
	for(int i=3;i>0;i--) tmp = tmp.squ();
	return _mul( tmp , x_2_3_2 );
}
};



template <typename gf_t>
struct fermat_inv<gf_t, 127> {
static const gf_t exec( const gf_t & x ) {
	/// Extend X to X^( 2^127 - 2 )
	gf_t x2 = x.squ();
	gf_t x3 = _mul( x2 , x );

	gf_t tmp = x3;
	tmp = tmp.squ();
	gf_t x_2_3_2 = tmp; ///
	tmp = tmp.squ();
	gf_t x_2_4_1 = _mul( tmp,x3); /// F

	tmp = x_2_4_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_8_1 = _mul( tmp,x_2_4_1); /// FF

	tmp = x_2_8_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_16_1 = _mul( tmp,x_2_8_1);  /// FFFF

	tmp = x_2_16_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_32_1 = _mul( tmp,x_2_16_1);  /// FFFF,FFFF

	tmp = x_2_32_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_64_1 = _mul( tmp,x_2_32_1);  /// FFFF,FFFF,FFFF,FFFF

	tmp = x_2_64_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_96_1 = _mul( tmp , x_2_32_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FFFF

	tmp = x_2_96_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_112_1 = _mul( tmp , x_2_16_1); /// FFFFx6,FFFF

	tmp = x_2_112_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_120_1 = _mul( tmp , x_2_8_1); /// FFFFx6,FFFF,FF

	tmp = x_2_120_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_124_1 = _mul( tmp , x_2_4_1); /// FFFFx6,FFFF,FFF

	tmp = x_2_124_1;
	for(int i=3;i>0;i--) tmp = tmp.squ();
	return _mul( tmp , x_2_3_2 );
}
};


template <typename gf_t>
struct fermat_inv<gf_t, 128> {
static const gf_t exec( const gf_t & x ) {
	/// Extend X to X^( 2^128 - 2 )
	gf_t x2 = x.squ();
	gf_t x3 = _mul( x2 , x );

	gf_t tmp = x3;
	for(int i=2;i>0;i--) tmp = tmp.squ();
	gf_t x_2_4_1 = _mul( tmp,x3); /// F
	gf_t x_2_4_2 = _mul( tmp,x2); /// E

	tmp = x_2_4_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_8_1 = _mul( tmp,x_2_4_1); /// FF

	tmp = x_2_8_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_16_1 = _mul( tmp,x_2_8_1);  /// FFFF

	tmp = x_2_16_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_32_1 = _mul( tmp,x_2_16_1);  /// FFFF,FFFF

	tmp = x_2_32_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_64_1 = _mul( tmp,x_2_32_1);  /// FFFF,FFFF,FFFF,FFFF

	tmp = x_2_64_1;
	for(int i=8*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_96_1 = _mul( tmp , x_2_32_1); /// FFFF,FFFF,FFFF,FFFF,FFFF,FFFF

	tmp = x_2_96_1;
	for(int i=4*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_112_1 = _mul( tmp , x_2_16_1); /// FFFFx6,FFFF

	tmp = x_2_112_1;
	for(int i=2*4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_120_1 = _mul( tmp , x_2_8_1); /// FFFFx6,FFFF,FF

	tmp = x_2_120_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	gf_t x_2_124_1 = _mul( tmp , x_2_4_1); /// FFFFx6,FFFF,FFF

	tmp = x_2_124_1;
	for(int i=4;i>0;i--) tmp = tmp.squ();
	return _mul( tmp , x_2_4_2 );
}
};


/*
template <> inline
gf2ext_u64<96> _inv( const gf2ext_u64<96> & x )
{
	static const unsigned W = 96;
	typedef gf2ext_u64<W> gf_t;
#ifdef CONFIG_PROFILE
count_inv();
#endif
	return fermat_inv<gf_t,96>::exec(x);
}
*/




#include "run_config.h"


#if defined( CONFIG_HAS_PCLMULQDQ ) || defined( CONFIG_PSHUFB )
#include "gf2ext-sse.hpp"
#endif

#if defined( CONFIG_NEON )
#include "gf2ext-neon.hpp"
#endif

#if defined( CONFIG_HAS_PCLMULQDQ )
#define GF2EXT gf2ext_sse
#elif defined( CONFIG_PSHUFB )
#define GF2EXT gf2ext_sse_tbl
#elif defined( CONFIG_NEON )
#define GF2EXT gf2ext_neon
#else
#define GF2EXT gf2ext_u64
#endif



#endif /// _GF2EXT_HPP_

