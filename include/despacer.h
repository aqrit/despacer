size_t despace_simple( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	for( ; length != 0; length-- ){
		uint8_t c = *src++;
		if( (c != 0x20) && (c != 0x0A) && (c != 0x0D) && (c != 0x09) ){
			*dst++ = c;
		}
	}
	return (size_t)(dst - ((uint8_t*)dst_void));
}


// make the store unconditional
size_t despace_branchless( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	for( ; length != 0; length-- ){
		uint8_t c = *src++;
		*dst = c;
		dst += !!((c != 0x20) && (c != 0x0A) && (c != 0x0D) && (c != 0x09));
	}
	return (size_t)(dst - ((uint8_t*)dst_void));
}


// explicit branchless
size_t despace_cmov( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	for( ; length != 0; length-- ){
		uint8_t c = *src++;
		uint64_t m = ( c > 0x20 ) ? 0xFFFFFFFFFFFFFFFF : 0xFFFFFFFEFFFFD9FF;
		*dst = c;
		dst += ((m >> (c & 63)) & 1);
	}
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_table( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	static const uint8_t table[256] = { 
		1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 
	};
	for( ; length != 0; length-- ){
		size_t c = *src++;
		*dst = (uint8_t)c;
		dst += table[c];
	}
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_block_simple( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	const uint64_t mask_09 = 0x0909090909090909;
	const uint64_t mask_0A = 0x0A0A0A0A0A0A0A0A;
	const uint64_t mask_0D = 0x0D0D0D0D0D0D0D0D;
	const uint64_t mask_20 = 0x2020202020202020;
	const uint64_t mask_7F = 0x7F7F7F7F7F7F7F7F;
	
	for( ; length >= 8; length-=8 ){
		uint64_t asrc = *((uint64_t*)src);
		
		uint64_t mask = asrc & mask_7F;
		mask =
			((mask ^ mask_09) + mask_7F) &
			((mask ^ mask_20) + mask_7F) &
			((mask ^ mask_0A) + mask_7F) &
			((mask ^ mask_0D) + mask_7F);
		mask = (mask | asrc) >> 7;
		// mask = bit0 of each byte is set if non-space

		for( int i = 0; i < 8; i++ ){
			*dst = src[i];
			dst += mask & 1;
			mask >>= 8;
		}
		src += 8;
	}

	// do remaining bytes
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


// for sparse whitespace only...
size_t despace_block_mux( void* dst_void, void* src_void, size_t length )
{
	uint64_t* src = (uint64_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const uint64_t mask_09 = 0x0909090909090909;
	const uint64_t mask_0A = 0x0A0A0A0A0A0A0A0A;
	const uint64_t mask_0D = 0x0D0D0D0D0D0D0D0D;
	const uint64_t mask_20 = 0x2020202020202020;
	const uint64_t mask_7F = 0x7F7F7F7F7F7F7F7F;

	for( ; length >= 8; length-=8 ){
		uint64_t asrc = *src++;
		size_t ws_cnt = 0;

		uint64_t mask = asrc & mask_7F;
		mask =
			((mask ^ mask_09) + mask_7F) &
			((mask ^ mask_20) + mask_7F) &
			((mask ^ mask_0A) + mask_7F) &
			((mask ^ mask_0D) + mask_7F);
		mask = ~mask & ~asrc & ~mask_7F;
		// mask = bit7 of each byte is set if space

		if(mask != 0){
			do{
				uint64_t pattern = (mask ^ (mask - 1)); // lowest set bit and below
				asrc = (((asrc << 8) & pattern) | (asrc & ~pattern)); // blend
				ws_cnt++; // count spaces
				mask &= ~pattern; // zero lowest set bit in mask
			}while(mask != 0);
			asrc >>= ws_cnt*8; // little endian
		}

		*((uint64_t*)dst) = asrc;
		dst += (8 - ws_cnt);
	}
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_sse2_detect( void* dst_void, void* src_void, size_t length )
{
	uint8_t* dst = (uint8_t*)dst_void;
	uint8_t* src = (uint8_t*)src_void;
	const __m128i mask_09 = _mm_set1_epi8(0x09);
	const __m128i mask_0A = _mm_set1_epi8(0x0A);
	const __m128i mask_0D = _mm_set1_epi8(0x0D);
	const __m128i mask_20 = _mm_set1_epi8(0x20);

	for( ; length >= 16; length-=16 ){
		__m128i v = _mm_loadu_si128((__m128i*)src);
		int m = _mm_movemask_epi8(
			_mm_or_si128(_mm_cmpeq_epi8(v, mask_09),
			_mm_or_si128(_mm_cmpeq_epi8(v, mask_0A),
			_mm_or_si128(_mm_cmpeq_epi8(v, mask_0D),
			_mm_cmpeq_epi8(v, mask_20)))));

		m = ~m;
		for( int i = 0; i < 16; i++ ){
			*dst = *src++;
			dst += m & 1;
			m >>= 1;
		}
	}
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_ssse3_cumsum( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m128i is_3or7 = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00);
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x01, 0x01, 0x00, 0x00, 0x01, 0x00, 0x00);
	const __m128i id = _mm_setr_epi8(
		0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
		0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F);
	const __m128i mask_02 = _mm_set1_epi8(0x02);
	const __m128i mask_04 = _mm_add_epi8(mask_02, mask_02);
	const __m128i mask_20 = _mm_slli_epi64(mask_02, 4);
	const __m128i mask_70 = _mm_set1_epi8(0x70);

	for( uint8_t* end = &src[(length & ~15)]; src != end; src += 16){
		__m128i a,b,c,d,s,t,v;
		size_t cnt0, cnt1;

		// load
		v = _mm_loadu_si128((__m128i*)src);

		// detect spaces ( 0x01 == space, 0x00 == non-space )
		s = _mm_or_si128(_mm_abs_epi8(_mm_cmpeq_epi8(mask_20, v)),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v)));

		// create non-space mask ( 0x00 == space, 0xFF == non-space )
		b = _mm_cmpeq_epi8(_mm_setzero_si128(), s);

		// (qword) prefix sum of spaces
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 8));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 16));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 32));

		// get non-space byte totals
		t = _mm_srli_epi64(s, 56); // hi-byte is total_spaces
		cnt0 = (uint32_t)_mm_cvtsi128_si32(t);
		t = _mm_unpackhi_epi64(t, t);
		cnt1 = (uint32_t)_mm_cvtsi128_si32(t);

		// compress
		b = _mm_andnot_si128(b, s); // zero non-spaces
		//
		c = _mm_srli_epi64(_mm_and_si128(mask_02, b), 9);
		d = _mm_srli_epi64(_mm_shuffle_epi8(is_3or7, b), 16);
		a = _mm_or_si128(_mm_cmpgt_epi8(mask_04, s), _mm_cmpeq_epi8(b, mask_04)); // match first 4 and below
		//
		s = _mm_add_epi8(s, c);
		s = _mm_add_epi8(s, d);
		//
		s = _mm_max_epu8(s, _mm_srli_epi64(_mm_andnot_si128(a, s), 32));
		v = _mm_shuffle_epi8(v, _mm_add_epi8(s, id));

		// store
		_mm_storel_epi64((__m128i*)dst, v);
		dst += 8 - cnt0;
		_mm_storel_epi64((__m128i*)dst, _mm_unpackhi_epi64(v, v));
		dst += 8 - cnt1;
	}
	dst += despace_branchless(dst, src, length & 15);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_ssse3_lut_1kb( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	static const uint64_t table[128] __attribute__((aligned(64))) = {
		0x0706050403020100, 0x0007060504030201, 0x0107060504030200, 0x0100070605040302,
		0x0207060504030100, 0x0200070605040301, 0x0201070605040300, 0x0201000706050403,
		0x0307060504020100, 0x0300070605040201, 0x0301070605040200, 0x0301000706050402,
		0x0302070605040100, 0x0302000706050401, 0x0302010706050400, 0x0302010007060504,
		0x0407060503020100, 0x0400070605030201, 0x0401070605030200, 0x0401000706050302,
		0x0402070605030100, 0x0402000706050301, 0x0402010706050300, 0x0402010007060503,
		0x0403070605020100, 0x0403000706050201, 0x0403010706050200, 0x0403010007060502,
		0x0403020706050100, 0x0403020007060501, 0x0403020107060500, 0x0403020100070605,
		0x0507060403020100, 0x0500070604030201, 0x0501070604030200, 0x0501000706040302,
		0x0502070604030100, 0x0502000706040301, 0x0502010706040300, 0x0502010007060403,
		0x0503070604020100, 0x0503000706040201, 0x0503010706040200, 0x0503010007060402,
		0x0503020706040100, 0x0503020007060401, 0x0503020107060400, 0x0503020100070604,
		0x0504070603020100, 0x0504000706030201, 0x0504010706030200, 0x0504010007060302,
		0x0504020706030100, 0x0504020007060301, 0x0504020107060300, 0x0504020100070603,
		0x0504030706020100, 0x0504030007060201, 0x0504030107060200, 0x0504030100070602,
		0x0504030207060100, 0x0504030200070601, 0x0504030201070600, 0x0504030201000706,
		0x0607050403020100, 0x0600070504030201, 0x0601070504030200, 0x0601000705040302,
		0x0602070504030100, 0x0602000705040301, 0x0602010705040300, 0x0602010007050403,
		0x0603070504020100, 0x0603000705040201, 0x0603010705040200, 0x0603010007050402,
		0x0603020705040100, 0x0603020007050401, 0x0603020107050400, 0x0603020100070504,
		0x0604070503020100, 0x0604000705030201, 0x0604010705030200, 0x0604010007050302,
		0x0604020705030100, 0x0604020007050301, 0x0604020107050300, 0x0604020100070503,
		0x0604030705020100, 0x0604030007050201, 0x0604030107050200, 0x0604030100070502,
		0x0604030207050100, 0x0604030200070501, 0x0604030201070500, 0x0604030201000705,
		0x0605070403020100, 0x0605000704030201, 0x0605010704030200, 0x0605010007040302,
		0x0605020704030100, 0x0605020007040301, 0x0605020107040300, 0x0605020100070403,
		0x0605030704020100, 0x0605030007040201, 0x0605030107040200, 0x0605030100070402,
		0x0605030207040100, 0x0605030200070401, 0x0605030201070400, 0x0605030201000704,
		0x0605040703020100, 0x0605040007030201, 0x0605040107030200, 0x0605040100070302,
		0x0605040207030100, 0x0605040200070301, 0x0605040201070300, 0x0605040201000703,
		0x0605040307020100, 0x0605040300070201, 0x0605040301070200, 0x0605040301000702,
		0x0605040302070100, 0x0605040302000701, 0x0605040302010700, 0x0605040302010007};

	const __m128i mask_01 = _mm_abs_epi8(_mm_cmpeq_epi8(_mm_setzero_si128(),_mm_setzero_si128()));
	const __m128i mask_20 = _mm_slli_epi64(mask_01, 5);
	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);

	for( uint8_t* end = &src[(length & ~15)]; src != end; src += 16){
		__m128i vector0 = _mm_loadu_si128((__m128i*)src);
		__m128i vector1 = _mm_shuffle_epi32( vector0, 0x0E );

		__m128i bytemask0 = _mm_or_si128(_mm_cmpeq_epi8(mask_20, vector0),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, vector0)));

		uint32_t bitmask0 = _mm_movemask_epi8(bytemask0) & 0x7F7F;
		__m128i hsum = _mm_sad_epu8(_mm_add_epi8(bytemask0, mask_01), _mm_setzero_si128());

		vector0 = _mm_shuffle_epi8(vector0, _mm_loadl_epi64((__m128i*) &table[(uint8_t)bitmask0]));
		_mm_storel_epi64((__m128i*)dst, vector0);
		dst += (uint32_t)_mm_cvtsi128_si32(hsum);

		vector1 = _mm_shuffle_epi8(vector1, _mm_loadl_epi64((__m128i*) &table[bitmask0 >> 8]));
		_mm_storel_epi64((__m128i*)dst, vector1);
		dst += (uint32_t)_mm_cvtsi128_si32(_mm_unpackhi_epi64(hsum, hsum));
	}

	dst += despace_branchless(dst, src, length & 15);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


__m128i lut_1mb [0x10000];
void gen_table_1mb( void ){
	for( size_t i = 0; i < 0x10000; i++ ){
		uint8_t* p = (uint8_t*) &lut_1mb[i];
		for( size_t j = 0; j < 16; j++ ){
			if( !(i & (1 << j)) ){
				*p++ = j;
			}
		}
	}
}
size_t despace_ssse3_lut_1mb( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m128i mask_01 = _mm_abs_epi8(_mm_cmpeq_epi8(_mm_setzero_si128(),_mm_setzero_si128()));
	const __m128i mask_20 = _mm_slli_epi64(mask_01, 5);
	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);

	for( uint8_t* end = &src[(length & ~15)]; src != end; src += 16){
		__m128i v = _mm_loadu_si128((__m128i*)src);
		__m128i bytemask = _mm_or_si128(_mm_cmpeq_epi8(mask_20, v),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v)));
		uint32_t bitmask = _mm_movemask_epi8(bytemask);
		v = _mm_shuffle_epi8(v, _mm_load_si128(&lut_1mb[bitmask]));
		_mm_storeu_si128((__m128i*)dst, v);
		__m128i hsum = _mm_sad_epu8(_mm_add_epi8(bytemask, mask_01), _mm_setzero_si128());
		dst += ((uint32_t)_mm_cvtsi128_si32(hsum)) + ((uint32_t)_mm_cvtsi128_si32(_mm_unpackhi_epi64(hsum, hsum)));
	}
	dst += despace_branchless(dst, src, length & 15);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_avx2_vpermd( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m256i mask_20  = _mm256_set1_epi8( 0x20 );
	const __m256i mask_70  = _mm256_set1_epi8( 0x70 );
	const __m256i lut_cntrl = _mm256_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00,
		//
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00
	);

	const __m256i mask_index = _mm256_set1_epi64x( 0x80A0908884828180 );
	const __m256i mask_shift = _mm256_set1_epi64x( 0x0000000008080808 );
	const __m256i mask_invert = _mm256_set1_epi64x( 0x0020100800000000 );
	const __m256i mask_fixup = _mm256_set_epi32(
		0x08080808, 0x0F0F0F0F, 0x00000000, 0x07070707,
		0x08080808, 0x0F0F0F0F, 0x00000000, 0x07070707
	);
	const __m256i lut = _mm256_set_epi32(
		0x04050607, // 0x03020100', 0x000000'07
		0x04050704, // 0x030200'00, 0x0000'0704
		0x04060705, // 0x030100'00, 0x0000'0705
		0x04070504, // 0x0300'0000, 0x00'070504
		0x05060706, // 0x020100'00, 0x0000'0706
		0x05070604, // 0x0200'0000, 0x00'070604
		0x06070605, // 0x0100'0000, 0x00'070605
		0x07060504  // 0x00'000000, 0x'07060504
	);

	for( uint8_t* end = &src[(length & ~31)]; src != end; src += 32){
		__m256i r0,r1,r2,r3,r4;

		r0 = _mm256_loadu_si256((__m256i *)src); // asrc

		r1 = _mm256_adds_epu8(mask_70, r0);
		r2 = _mm256_cmpeq_epi8(mask_20, r0);
		r1 = _mm256_shuffle_epi8(lut_cntrl, r1);
		r1 = _mm256_or_si256(r1, r2); // bytemask of spaces

		r2 = _mm256_andnot_si256(r1, mask_index);
		r1 = _mm256_and_si256(r1, mask_shift);
		r2 = _mm256_sad_epu8(r2, mask_invert); // bitmap[0:5], popcnt[7:15]
		r1 = _mm256_sad_epu8(r1, _mm256_setzero_si256()); // shift amount
		r3 = _mm256_slli_epi64(r2, 29); // move hi index to 2nd dword
		r4 = _mm256_srli_epi64(r2, 7); // popcnt
		r2 = _mm256_or_si256(r2, r3);
		r2 = _mm256_permutevar8x32_epi32(lut, r2);
		r2 = _mm256_xor_si256(r2, mask_fixup);
		r2 = _mm256_srlv_epi64(r2, r1);
		r0 = _mm256_shuffle_epi8(r0, r2);

		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 0);
		dst += _mm256_extract_epi64(r4, 0);
		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 1);
		dst += _mm256_extract_epi64(r4, 1);
		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 2);
		dst += _mm256_extract_epi64(r4, 2);
		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 3);
		dst += _mm256_extract_epi64(r4, 3);
	}
	dst += despace_branchless(dst, src, length & 31);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_avx2_lut_1mb( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m256i mask_70 = _mm256_set1_epi8( 0x70 );
	const __m256i mask_20 = _mm256_set1_epi8( 0x20 );
	const __m256i lut_cntrl = _mm256_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00,
		//
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00
		);

	for( uint8_t* end = &src[(length & ~31)]; src != end; src += 32){
		__m256i v = _mm256_loadu_si256((__m256i*)src);
		__m256i bytemask = _mm256_or_si256(_mm256_cmpeq_epi8(mask_20, v),
			_mm256_shuffle_epi8(lut_cntrl, _mm256_adds_epu8(mask_70, v)));
		uint32_t bitmask = _mm256_movemask_epi8(bytemask);
		_mm_storeu_si128((__m128i*)dst, _mm_shuffle_epi8(_mm256_castsi256_si128(v), _mm_load_si128(&lut_1mb[bitmask & 0xFFFF])));
		_mm_storeu_si128((__m128i*)&dst[16 - _mm_popcnt_u32(bitmask & 0xFFFF)], _mm_shuffle_epi8(_mm256_extracti128_si256(v,1), _mm_load_si128(&lut_1mb[(bitmask >> 16)])));
		dst += (32 - _mm_popcnt_u32(bitmask));
	}
	dst += despace_branchless(dst, src, length & 31);
	return (size_t)(dst - ((uint8_t*)dst_void));
}
