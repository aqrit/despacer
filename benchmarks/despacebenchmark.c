#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <emmintrin.h> // sse2
#include <tmmintrin.h> // ssse3
#include <nmmintrin.h> // sse42

#include "despacer.h"


// generate random number from range 
// (rejection sampling to avoid modulo bias)
static uint32_t rand_range(uint32_t min, uint32_t max)
{
	assert(min <= max);
	assert((max-min) <= RAND_MAX);

	uint32_t r;
	uint32_t n = max - min;
	uint32_t m = n;
	m |= m >> 16; m |= m >> 8; m |= m >> 4; m |= m >> 2; m |= m >> 1;
	do{
		r = ((uint32_t)rand()) & m;
	}while(r > n);
	return r + min;
}


// randomly shuffle byte array (Knuth-Fisher–Yates shuffle)
static void rand_shuffle(uint8_t* arr, size_t n)
{
	assert(arr != 0);
	assert(n != 0);

	size_t i, j, tmp;
	for(i = n - 1; i != 0; i--){
		j = rand_range(0,i);
		tmp = arr[j];
		arr[j] = arr[i];
		arr[i] = tmp;
	}
}


// fill buf with random bytes and passed amount of whitespace
static void fillwithtext(uint8_t* buf, size_t size, size_t num_spaces)
{
	assert(buf != 0);
	assert(num_spaces <= size);

	size_t i;
	for(i = 0; i < num_spaces; i++){
		buf[i] = (uint8_t)(0x090A0D20 >> (rand() & 0x18));
	}
	for(; i < size; i++){
		uint8_t c;
		do{
			c = (uint8_t)rand();
		}while(c == 0x09 || c == 0x0A || c == 0x0D || c == 0x20);
		buf[i] = c; 
	}
	rand_shuffle(buf, size);	
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
	printf("%6I64d\n", t);
}


#define BUF_SIZE 1024
int main(int argc, char ** argv)
{
	(void)argc;
	(void)argv;

	srand(time(0));

	static uint8_t src[BUF_SIZE];
	static uint8_t dst[BUF_SIZE];
	static uint8_t ans[BUF_SIZE];
	size_t num_spaces = 128;
	assert(num_spaces <= BUF_SIZE);

	fillwithtext(src, BUF_SIZE, num_spaces);
	despace_simple(ans, src, BUF_SIZE);

	test_time( "despace_simple", &despace_simple,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_setcc", &despace_setcc,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_branchless", &despace_branchless,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_cmov", &despace_cmov,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_table", &despace_table,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_block_simple", &despace_block_simple,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_block_branchless", &despace_block_branchless,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_block_mux", &despace_block_mux,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_sse2_cumsum", &despace_sse2_cumsum,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_sse2_detect", &despace_sse2_detect,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_lut", &despace_ssse3_lut,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );

	return 0;
}