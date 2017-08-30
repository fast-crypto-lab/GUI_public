#ifndef _CONFIG_H_
#define _CONFIG_H_

#include "scheme.h"

//#define QUARTZ1273
//#define QUARTZ103
//#define QUARTZ95
//#define QUARTZ94
//#define QUARTZ96

#ifndef QUARTZ96
#ifndef QUARTZ59556
#ifndef QUARTZ1273
#ifndef QUARTZ127
#ifndef QUARTZ94
#ifndef QUARTZ95
#ifndef QUARTZ103
#ifndef QUARTZ128
#define QUARTZ96
//#define QUARTZ127
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif



#ifdef QUARTZ96
#define CORE_SIZE 96
#define MAX_DEG 5
#define VINEGAR 6
#define MINUS 6
#define REPEAT 3
/// XXX:
#define PUBKEY_BYTES 63036
#define SECKEY_BYTES 3175

#endif // end of QUARTZ59556


#ifdef QUARTZ59556
#define CORE_SIZE 95
#define MAX_DEG 5
#define VINEGAR 6
#define MINUS 5
#define REPEAT 3
/// XXX:
#define PUBKEY_BYTES 61812
#define SECKEY_BYTES 3150

#endif // end of QUARTZ59556


#ifdef QUARTZ1273
#define CORE_SIZE 127
#define MAX_DEG 9
#define VINEGAR 6
#define MINUS 4
#define REPEAT 3

#define PUBKEY_BYTES 142576
#define SECKEY_BYTES 5350

#endif // end of QUARTZ1273


#ifdef QUARTZ127
#define CORE_SIZE 127
#define MAX_DEG 9
#define VINEGAR 6
#define MINUS 4
#define REPEAT 4

#define PUBKEY_BYTES 142576
#define SECKEY_BYTES 5350

#endif // end of QUARTZ127


#ifdef QUARTZ128
#define CORE_SIZE 128
#define MAX_DEG 9
#define VINEGAR 6
#define MINUS 4
#define REPEAT 4

#define PUBKEY_BYTES 144720
#define SECKEY_BYTES 6336

#endif // end of QUARTZ128



#ifdef QUARTZ103
#define CORE_SIZE 103
#define MAX_DEG 129
#define VINEGAR 4
#define MINUS 3
#define REPEAT 4

#define PUBKEY_BYTES 75514
#define SECKEY_BYTES 3774

#endif // end of QUARTZ103

#ifdef QUARTZ95
#define CORE_SIZE 95
#define MAX_DEG 9
#define VINEGAR 5
#define MINUS 5
#define REPEAT 3

#define PUBKEY_BYTES 60600
#define SECKEY_BYTES 3053

#endif // end of QUARTZ95


#ifdef QUARTZ94
#define CORE_SIZE 94
#define MAX_DEG 17
#define VINEGAR 4
#define MINUS 4
#define REPEAT 4

#define PUBKEY_BYTES 58212
#define SECKEY_BYTES 2943
#endif // end of QUARTZ94


#define N (CORE_SIZE+VINEGAR)
#define M (CORE_SIZE-MINUS)

#define SECMSG_BYTES (((M+(MINUS+VINEGAR)*REPEAT)+7)/8)


#define S_WIDTH N
#define T_WIDTH CORE_SIZE

#define PUBKEY_NUM_TERMS ((N)*(N+1)/2)

#define SIZE_BYTE_M ((M+7)>>3)
#define SIZE_BYTE_N ((N+7)>>3)
#define SIZE_BYTE_CORE ((CORE_SIZE+7)>>3)




#endif
