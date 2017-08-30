#ifndef _STAT_PROFILE_H_
#define _STAT_PROFILE_H_

#include "benchmark.h"
#include <stdint.h>

void stat_reset();
void stat_say();


extern benchmark prepare_poly;
extern benchmark do_poly;
extern benchmark ext_gcd_poly;


void count_mul();
void count_squ();
void count_inv();
void count_poly_add();

//extern unsigned phase;

void profile_start();
void profile_end();
void new_phase();

uint64_t num_mul(unsigned);
uint64_t num_squ(unsigned);
uint64_t num_inv(unsigned);
uint64_t num_poly(unsigned);
//extern uint64_t num_mul[3];
//extern uint64_t num_squ[3];
//extern uint64_t num_inv[3];
//extern uint64_t num_poly[3];

#endif
