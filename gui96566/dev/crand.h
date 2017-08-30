#ifndef _CRAND_H_
#define _CRAND_H_

#include "run_config.h"

#if defined(_NO_OPENSSL_)

#include <stdint.h>
#include <stdlib.h>

static inline
void RAND_bytes( unsigned char* buff, unsigned s )
{
	union{
	unsigned char v8[4];
	uint32_t v32;
	};
	for(unsigned i=0;i<s;i++){
		if( 0==(i&1) ) v32 = rand();
		buff[i]=v8[i&1];
	}
}

#else

#include <openssl/rand.h>

#endif /* _NO_OPENSSL_ */

#endif
