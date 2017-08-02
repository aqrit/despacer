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
#include <immintrin.h> // avx2

#include "despacer.h"
#include "useless.h"


static void fillwithtext(uint8_t* buf, size_t size)
{
	uint8_t space[] = { 0x09, 0x0A, 0x0D, 0x20 };
	uint8_t nonspace[] = { 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
		0x38, 0x39, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46 };

	for( size_t i = 0; i < size >> 4; i++ ){
		for( size_t j = 0; j < 16; j++ ){
			buf[i*16+j] = (i & (1 << j)) ?  space[rand() & 0x03] : nonspace[rand() & 0x0F];
		}
	}
	for( size_t i = 0; i < size % 16; i++ ) buf[(size >> 4) * 16 + i] = 0x46;

	if( size < 16 ) return;

	for(size_t i = (size >> 4) - 1; i != 0; i--){
		uint32_t j = 0;
		uint32_t m = i;
		m |= m >> 16; m |= m >> 8; m |= m >> 4; m |= m >> 2; m |= m >> 1;
		do{
			j = ((rand() << 15) ^ rand()) & m;
		}while(j > i);

		__m128i tmp_j = _mm_loadu_si128((__m128i*)&buf[j]);
		__m128i tmp_i = _mm_loadu_si128((__m128i*)&buf[i]);
		_mm_storeu_si128((__m128i*)&buf[j], tmp_i);
		_mm_storeu_si128((__m128i*)&buf[i], tmp_j);
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
	printf("%10" PRIu64 "\n", t);
}


#define BUF_SIZE 0x8000000
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
	printf( "buf_size: %d  num_spaces: %" PRIu64 "\n\n", BUF_SIZE, num_spaces );

	gen_table_1mb();

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
	test_time( "despace_ssse3_cumsum", &despace_ssse3_cumsum,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_lut_1kb", &despace_ssse3_lut_1kb,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_lut_1mb", &despace_ssse3_lut_1mb,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_avx2_lut_1mb", &despace_avx2_lut_1mb,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_avx2_vpermd", &despace_avx2_vpermd,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );

	printf("\n%24s\n","---useless---");
	gen_table_1mb_2();
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
	test_time( "despace_ssse3_lut_512b", &despace_ssse3_lut_512b,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );
	test_time( "despace_ssse3_lut_1mb_2", &despace_ssse3_lut_1mb_2,
		dst, src, BUF_SIZE, ans, BUF_SIZE - num_spaces );

	return 0;
}
