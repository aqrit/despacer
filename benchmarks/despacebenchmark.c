#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <emmintrin.h> // sse2
#include <tmmintrin.h> // ssse3
#include <smmintrin.h> // sse41
#include <nmmintrin.h> // sse42

#include "despacer.h"
#include "useless.h"


// fill buf with bytes... coinflip to choose either 0x41 or 0x20 
static void fillwithtext(uint8_t* buf, size_t size)
{
	int n = 0;
	uint64_t r = 0;
	for( size_t i = 0; i < size; i++ ){
		if( --n < 0 ){ r = rand(); n = 14; }
		buf[i] = (r & (1 << n)) ? 0x41 : 0x20;
	}
}


#define RDTSC_START(cycles)                                     \
  do{                                                           \
    register unsigned cyc_high, cyc_low;                        \
    __asm volatile("cpuid\n\t"                                  \
                   "rdtsc\n\t"                                  \
                   "cpuid\n\t"                                  \
                   "rdtsc\n\t"                                  \
                   "cpuid\n\t"                                  \
                   "rdtsc\n\t"                                  \
                   "mov %%edx, %0\n\t"                          \
                   "mov %%eax, %1\n\t"                          \
                   : "=r"(cyc_high), "=r"(cyc_low)              \
                   ::"%rax", "%rbx", "%rcx", "%rdx", "memory"); \
    (cycles) = ((uint64_t)cyc_high << 32) | cyc_low;            \
  }while(0)

#define RDTSC_FINAL(cycles)                                     \
  do{                                                           \
    register unsigned cyc_high, cyc_low;                        \
    __asm volatile("rdtscp\n\t"                                 \
                   "mov %%edx, %0\n\t"                          \
                   "mov %%eax, %1\n\t"                          \
                   "cpuid\n\t"                                  \
                   : "=r"(cyc_high), "=r"(cyc_low)              \
                   ::"%rax", "%rbx", "%rcx", "%rdx", "memory"); \
    (cycles) = ((uint64_t)cyc_high << 32) | cyc_low;            \
  }while(0)


// benchmark
static uint64_t best_time(
	size_t (*fn_ptr)(void*,void*,size_t),
	uint8_t* out,
	const uint8_t* const in,
	size_t length )
{
	uint64_t tm1, tm2;
	uint64_t min_diff = 0xFFFFFFFFFFFFFFFF;
	for(int i = 0; i < 100; i++){
		memcpy(out, in, length);
		RDTSC_START(tm1);
		(*fn_ptr)(out, out, length);
		RDTSC_FINAL(tm2);
		uint64_t tmus = tm2 - tm1;
		if(tmus < min_diff) min_diff = tmus;
	}
	return min_diff;
}


//
static void test_time(
	char* function_name,
	size_t (*fn_ptr)(void*,void*,size_t), 
	uint8_t* out,
	const uint8_t* const in,
	size_t in_len,
	const uint8_t* const ans,
	size_t ans_len )
{
	printf("%24s: ", function_name); fflush(0);
	memcpy(out, in, in_len);
	size_t n = (*fn_ptr)(out, out, in_len);
	if((n != ans_len) || (memcmp(out,ans,n) != 0)){
		printf("BUG!!! "); fflush(0);
	}
	uint64_t t = best_time(fn_ptr, out, in, in_len);
	printf("%5" PRIu64 "\n", t);
}


#define BUF_SIZE 0x1000000
int main(int argc, char ** argv)
{
	(void)argc;
	(void)argv;

	srand(time(0));

	static uint8_t src[BUF_SIZE];
	static uint8_t dst[BUF_SIZE];
	static uint8_t ans[BUF_SIZE];
	size_t num_spaces;
	assert(num_spaces <= BUF_SIZE);

	fillwithtext(src, BUF_SIZE);
	num_spaces = BUF_SIZE - despace_simple(ans, src, BUF_SIZE);
	gen_table_1mb();

	/*
	test_time( "despace_simple", &despace_simple,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_branchless", &despace_branchless,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_cmov", &despace_cmov,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_table", &despace_table,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_block_simple", &despace_block_simple,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_block_mux", &despace_block_mux,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_sse2_detect", &despace_sse2_detect,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	*/
	test_time( "despace_ssse3_cumsum", &despace_ssse3_cumsum,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_lut_512b", &despace_ssse3_lut_512b,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_lut_1kb", &despace_ssse3_lut_1kb,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_lut_1mb", &despace_ssse3_lut_1mb,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );

	/*
	printf("\n%24s\n","---useless---");
	test_time( "despace_setcc", &despace_setcc,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_block_shift", &despace_block_shift,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_sse42_detect", &despace_sse42_detect,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_sse42_scan", &despace_sse42_scan,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_detect", &despace_ssse3_detect,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_sse2_cumsum", &despace_sse2_cumsum,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_sse41_cumsum", &despace_sse41_cumsum,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	*/
	return 0;
}
