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

	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);
	const __m128i id = _mm_setr_epi8(
		0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
		0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F);
	const __m128i mask_02 = _mm_set1_epi8(0x02);
	const __m128i mask_04 = _mm_add_epi8(mask_02, mask_02);
	const __m128i mask_20 = _mm_slli_epi64(mask_02, 4);
	const __m128i mask_70 = _mm_set1_epi8(0x70);

	for( ; length >= 16; length-=16 ){
		// load
		__m128i v = _mm_loadu_si128((__m128i*)src);
		src += 16;

		// detect spaces
		__m128i ws_FF = _mm_or_si128(_mm_cmpeq_epi8(mask_20, v),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v)));

		// (qword) prefix sum of spaces
		__m128i ws_01 = _mm_abs_epi8( ws_FF );
		__m128i s = ws_01;
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 8));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 16));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 32));

		// get non-space byte totals
		__m128i cnt = _mm_srli_epi64(s, 56);
		size_t popcnt0 = 8 - _mm_cvtsi128_si32(cnt);
		size_t popcnt1 = 8 - _mm_cvtsi128_si32(_mm_unpackhi_epi64(cnt, cnt));

		//
		s = _mm_andnot_si128(_mm_add_epi8(ws_FF, _mm_and_si128(ws_01, s)), s);
		s = _mm_max_epu8(s, _mm_srli_epi64(_mm_and_si128(_mm_cmpgt_epi8(_mm_and_si128(mask_02, s), _mm_setzero_si128()), s), 16));
		s = _mm_max_epu8(s, _mm_srli_epi64(_mm_andnot_si128(_mm_cmpgt_epi8(mask_04, s), s), 32));
		v = _mm_shuffle_epi8(v, _mm_add_epi8(s, id));

		// store
		_mm_storel_epi64((__m128i*)dst, v);
		dst += popcnt0;
		_mm_storel_epi64((__m128i*)dst, _mm_unpackhi_epi64(v, v));
		dst += popcnt1;
	}
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_ssse3_lut_1kb( void* dst_void, void* src_void, size_t length )
{
	__m128i * src = (__m128i *)src_void;
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

	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i mask_20 = _mm_set1_epi8( 0x20 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);

	for( ; length >= 16; length-=16 ){
		// load
		__m128i vector0 = _mm_loadu_si128(src);
		src++;
		__m128i vector1 = _mm_unpackhi_epi64(vector0, vector0);

		// detect spaces
		__m128i bytemask0 = _mm_or_si128(_mm_cmpeq_epi8(mask_20, vector0),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, vector0)));

		// sort (compress)
		int bitmask0 = _mm_movemask_epi8(bytemask0);
		vector0 = _mm_shuffle_epi8(vector0, _mm_loadl_epi64((__m128i*) &table[bitmask0 & 0x7F]));
		vector1 = _mm_shuffle_epi8(vector1, _mm_loadl_epi64((__m128i*) &table[(bitmask0 >> 8) & 0x7F]));

		// count non-spaces
		__m128i hsum0 = _mm_sad_epu8( _mm_setzero_si128(), bytemask0 );
		__m128i hsum1 = _mm_unpackhi_epi64(hsum0, hsum0);
		size_t popcnt0 = ((uint8_t)( 8 + _mm_cvtsi128_si32(hsum0)));
		size_t popcnt1 = ((uint8_t)( 8 + _mm_cvtsi128_si32(hsum1)));

		// store
		_mm_storel_epi64((__m128i*)dst, vector0);
		dst += popcnt0;
		_mm_storel_epi64((__m128i*)dst, vector1);
		dst += popcnt1;
	}

	// do remaining bytes
	dst += despace_branchless(dst, src, length);

	return (size_t)(dst - ((uint8_t*)dst_void));
}
