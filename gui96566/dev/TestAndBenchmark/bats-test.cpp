

#include "quartz.hpp"

#include "quartz.h"


#define TEST_RUN 50

int main()
{
	printf("QUARTZ<%d,%d,%d,%d> %d\n",CORE_SIZE,MAX_DEG,MINUS,VINEGAR,REPEAT);

	printf("signature byte: %d\n", vec_sign_t::num_byte());

	printf("pub key: %d,\n", quartz_pub_key_t::num_byte(),PUBLICKEY_BYTES);

	printf("sec key: %d,\n", quartz_sec_key_t::num_byte(),SECRETKEY_BYTES);

	printf("sec meg: %d\n", SIGNATURE_BYTES );

	unsigned char m[32];
	unsigned char sm[SIGNATURE_BYTES];

	unsigned char pk[PUBLICKEY_BYTES];
	unsigned char sk[SECRETKEY_BYTES];
	unsigned long long pklen;
	unsigned long long sklen;
	unsigned long long smlen;
	unsigned long long mlen = 32;

	int kp = keypair( sk , &sklen , pk , &pklen );

	RAND_bytes( m , 32 );

	int ssh = signatureofshorthash(sm,&smlen,m,mlen,sk,sklen);
	printf("smlen: %d\n", smlen);

	int vf = verification(m,mlen,sm,smlen,pk,pklen);

	printf("verify: %d\n",vf);


	int succ = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		RAND_bytes( m , 32 );
		signatureofshorthash(sm,&smlen,m,mlen,sk,sklen);
		vf = verification(m,mlen,sm,smlen,pk,pklen);
		if( 0 == vf ) succ++;
	}
	printf("%d/%d passed.\n", succ , TEST_RUN);
	printf("\n\n");

	return 0;
}
