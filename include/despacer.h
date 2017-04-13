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


size_t despace_setcc( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	for( ; length != 0; length-- ){
		uint8_t c = *src++;
		*dst = c;
		dst += ((c == 0x20) | (c == 0x0A) | (c == 0x0D)) ^ (c != 0x09);
	}
	return (size_t)(dst - ((uint8_t*)dst_void));
}


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


// branchless explicit
size_t despace_cmov( void* dst_void, void* src_void, size_t length )
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	for( ; length != 0; length-- ){
		uint8_t c = *src++;
		uint64_t m = ( c > 0x20 ) ? 0xFFFFFFFFFFFFFFFF : 0xFFFFFFFEFFFFD9FF;
		*dst = c;
		dst += ((m >> (c & 63)) & 1); 
		// note: the bitwise-and is needed to avoid "undefined behavior"
		// unfortunaly, it isn't being optimized away...
		// but here there seems to be little impact on speed.
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
	uint64_t* src = (uint64_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	const uint64_t mask_09 = 0x0909090909090909;
	const uint64_t mask_0A = 0x0A0A0A0A0A0A0A0A;
	const uint64_t mask_0D = 0x0D0D0D0D0D0D0D0D;
	const uint64_t mask_20 = 0x2020202020202020;
	const uint64_t mask_7F = 0x7F7F7F7F7F7F7F7F;
	
	for( ; length >= 8; length-=8 ){
		uint64_t asrc = *src++;
		
		uint64_t mask = asrc & mask_7F;
		mask = 
			((mask ^ mask_09) + mask_7F) &
			((mask ^ mask_20) + mask_7F) &
			((mask ^ mask_0A) + mask_7F) &
			((mask ^ mask_0D) + mask_7F);
		mask = (( mask | asrc ) & ~mask_7F );
		// mask = bit7 of each byte is set if non-space

		for( ; mask != 0; mask >>= 8 ){
			if( ((uint8_t)mask) != 0 ){
				*dst = (uint8_t)asrc;
				dst++;
			}
			asrc >>= 8;
		}
	}

	// do remaining bytes
	dst += despace_branchless(dst, src, length);

	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_block_branchless( void* dst_void, void* src_void, size_t length )
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

		uint64_t mask = asrc & mask_7F;
		mask = 
			((mask ^ mask_09) + mask_7F) &
			((mask ^ mask_20) + mask_7F) &
			((mask ^ mask_0A) + mask_7F) &
			((mask ^ mask_0D) + mask_7F);
		mask = (( mask | asrc ) & ~mask_7F ) >> 7;
		// mask = bit0 of each byte is set if non-space

		*dst = ((uint8_t)asrc);
		dst += ((uint8_t)mask);
		*dst = ((uint8_t)(asrc >> 8));
		dst += ((uint8_t)(mask >> 8));
		asrc >>= 16;
		mask >>= 16;
		*dst = ((uint8_t)asrc);
		dst += ((uint8_t)mask);
		*dst = ((uint8_t)(asrc >> 8));
		dst += ((uint8_t)(mask >> 8));
		asrc >>= 16;
		mask >>= 16;
		*dst = ((uint8_t)asrc);
		dst += ((uint8_t)mask);
		*dst = ((uint8_t)(asrc >> 8));
		dst += ((uint8_t)(mask >> 8));
		asrc >>= 16;
		mask >>= 16;
		*dst = ((uint8_t)asrc);
		dst += ((uint8_t)mask);
		*dst = ((uint8_t)(asrc >> 8));
		dst += ((uint8_t)(mask >> 8));
	}

	// do remaining bytes
	dst += despace_branchless(dst, src, length);

	return (size_t)(dst - ((uint8_t*)dst_void));
}


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

	// do remaining bytes
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

	// do remaining bytes
	dst += despace_branchless(dst, src, length);

	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_sse2_cumsum( void* dst_void, void* src_void, size_t length )
{
	__m128i* src = (__m128i*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m128i mask_09 = _mm_set1_epi8(0x09);
	const __m128i mask_0A = _mm_set1_epi8(0x0A);
	const __m128i mask_0D = _mm_set1_epi8(0x0D);
	const __m128i mask_20 = _mm_set1_epi8(0x20);
	const __m128i mask_01 = _mm_set1_epi8(0x01);

	for( ; length >= 16; length-=16 ){
		// load
		__m128i v = _mm_loadu_si128(src);
		src++;

		// detect spaces
		__m128i ws =
			_mm_or_si128(_mm_cmpeq_epi8(v, mask_09),
			_mm_or_si128(_mm_cmpeq_epi8(v, mask_0A),
			_mm_or_si128(_mm_cmpeq_epi8(v, mask_0D),
			_mm_cmpeq_epi8(v, mask_20))));

		// (qword) prefix sum of spaces
		__m128i s = _mm_and_si128(mask_01, ws);
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 8));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 16));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 32));

		//
		__m128i cnt = _mm_srli_epi64(s, 56);
		size_t popcnt0 = 8 - _mm_cvtsi128_si32(cnt);
		size_t popcnt1 = 8 - _mm_cvtsi128_si32(_mm_unpackhi_epi64(cnt, cnt));

		// zero space bytes in prefix sum
		s = _mm_andnot_si128(ws, s);

		//
		__m128i p,m,t;

		p = mask_01;
		t = _mm_srli_epi64(s, 8);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_or_si128(_mm_andnot_si128(m, s), _mm_and_si128(t, m));
		v = _mm_or_si128(_mm_andnot_si128(m, v), _mm_and_si128(_mm_srli_epi64(v, 8), m));

		p = _mm_add_epi8(p, p); // mask_02
		t = _mm_srli_epi64(s, 16);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_or_si128(_mm_andnot_si128(m, s), _mm_and_si128(t, m));
		v = _mm_or_si128(_mm_andnot_si128(m, v), _mm_and_si128(_mm_srli_epi64(v, 16), m));

		p = _mm_add_epi8(p, p); // mask_04
		t = _mm_srli_epi64(s, 32);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		// s = not needed on last pass
		v = _mm_or_si128(_mm_andnot_si128(m, v), _mm_and_si128(_mm_srli_epi64(v, 32), m));

		// store
		_mm_storel_epi64((__m128i*)dst, v);
		dst += popcnt0;
		_mm_storel_epi64((__m128i*)dst, _mm_unpackhi_epi64(v, v));
		dst += popcnt1;
	}

	// do remaining bytes
	dst += despace_branchless(dst, src, length);

	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_sse3_cumsum( void* dst_void, void* src_void, size_t length )
{
	__m128i* src = (__m128i*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i mask_20 = _mm_set1_epi8( 0x20 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);
	const __m128i id = _mm_setr_epi8(
		0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
		0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F);
	const __m128i mask_01 = _mm_set1_epi8(0x01);

	for( ; length >= 16; length-=16 ){
		// load
		__m128i v = _mm_loadu_si128(src);
		src++;

		// detect spaces
		__m128i ws = _mm_or_si128(_mm_cmpeq_epi8(mask_20, v),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v)));

		// (qword) prefix sum of spaces
		__m128i s = _mm_and_si128(mask_01, ws);
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 8));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 16));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 32));

		//
		__m128i cnt = _mm_srli_epi64(s, 56);
		size_t popcnt0 = 8 - _mm_cvtsi128_si32(cnt);
		size_t popcnt1 = 8 - _mm_cvtsi128_si32(_mm_unpackhi_epi64(cnt, cnt));

		// zero space bytes in prefix sum
		s = _mm_andnot_si128(ws, s);

		//
		__m128i p,m,t;

		p = mask_01;
		t = _mm_srli_epi64(s, 8);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_or_si128(_mm_andnot_si128(m, s), _mm_and_si128(t, m));

		p = _mm_add_epi8(p, p); // mask_02
		t = _mm_srli_epi64(s, 16);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_or_si128(_mm_andnot_si128(m, s), _mm_and_si128(t, m));

		p = _mm_add_epi8(p, p); // mask_04
		t = _mm_srli_epi64(s, 32);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_or_si128(_mm_andnot_si128(m, s), _mm_and_si128(t, m));

		v = _mm_shuffle_epi8(v, _mm_add_epi8(id, s));

		// store
		_mm_storel_epi64((__m128i*)dst, v);
		dst += popcnt0;
		_mm_storel_epi64((__m128i*)dst, _mm_unpackhi_epi64(v, v));
		dst += popcnt1;
	}

	// do remaining bytes
	dst += despace_branchless(dst, src, length);

	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_sse41_cumsum( void* dst_void, void* src_void, size_t length )
{
	__m128i* src = (__m128i*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i mask_20 = _mm_set1_epi8( 0x20 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);
	const __m128i id = _mm_setr_epi8(
		0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
		0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F);
	const __m128i mask_01 = _mm_set1_epi8(0x01);

	for( ; length >= 16; length-=16 ){
		// load
		__m128i v = _mm_loadu_si128(src);
		src++;

		// detect spaces
		__m128i ws = _mm_or_si128(_mm_cmpeq_epi8(mask_20, v),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v)));

		// (qword) prefix sum of spaces
		__m128i s = _mm_and_si128(mask_01, ws); 
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 8));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 16));
		s = _mm_add_epi8(s, _mm_slli_epi64(s, 32));

		//
		__m128i cnt = _mm_srli_epi64(s, 56);
		size_t popcnt0 = 8 - _mm_cvtsi128_si32(cnt);
		size_t popcnt1 = 8 - _mm_cvtsi128_si32(_mm_unpackhi_epi64(cnt, cnt));

		// zero space bytes in prefix sum
		s = _mm_andnot_si128(ws, s);

		//
		__m128i p,m,t;

		p = mask_01;
		t = _mm_srli_epi64(s, 8);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_blendv_epi8(s, t, m);

		p = _mm_add_epi8(p, p); // mask_02
		t = _mm_srli_epi64(s, 16);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_blendv_epi8(s, t, m);

		p = _mm_add_epi8(p, p); // mask_04
		t = _mm_srli_epi64(s, 32);
		m = _mm_cmpeq_epi8(p, _mm_and_si128(p, t));
		s = _mm_blendv_epi8(s, t, m);

		v = _mm_shuffle_epi8(v, _mm_add_epi8(id, s));

		// store
		_mm_storel_epi64((__m128i*)dst, v);
		dst += popcnt0;
		_mm_storel_epi64((__m128i*)dst, _mm_unpackhi_epi64(v, v));
		dst += popcnt1;
	}

	// do remaining bytes
	dst += despace_branchless(dst, src, length);

	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_ssse3_lut( void* dst_void, void* src_void, size_t length )
{
	__m128i * src = (__m128i *)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	static const uint8_t table[128][8] = {
		{0,1,2,3,4,5,6,7}, {1,2,3,4,5,6,7,0}, {0,2,3,4,5,6,7,1}, {2,3,4,5,6,7,0,1},
		{0,1,3,4,5,6,7,2}, {1,3,4,5,6,7,0,2}, {0,3,4,5,6,7,1,2}, {3,4,5,6,7,0,1,2},
		{0,1,2,4,5,6,7,3}, {1,2,4,5,6,7,0,3}, {0,2,4,5,6,7,1,3}, {2,4,5,6,7,0,1,3},
		{0,1,4,5,6,7,2,3}, {1,4,5,6,7,0,2,3}, {0,4,5,6,7,1,2,3}, {4,5,6,7,0,1,2,3},
		{0,1,2,3,5,6,7,4}, {1,2,3,5,6,7,0,4}, {0,2,3,5,6,7,1,4}, {2,3,5,6,7,0,1,4},
		{0,1,3,5,6,7,2,4}, {1,3,5,6,7,0,2,4}, {0,3,5,6,7,1,2,4}, {3,5,6,7,0,1,2,4},
		{0,1,2,5,6,7,3,4}, {1,2,5,6,7,0,3,4}, {0,2,5,6,7,1,3,4}, {2,5,6,7,0,1,3,4},
		{0,1,5,6,7,2,3,4}, {1,5,6,7,0,2,3,4}, {0,5,6,7,1,2,3,4}, {5,6,7,0,1,2,3,4},
		{0,1,2,3,4,6,7,5}, {1,2,3,4,6,7,0,5}, {0,2,3,4,6,7,1,5}, {2,3,4,6,7,0,1,5},
		{0,1,3,4,6,7,2,5}, {1,3,4,6,7,0,2,5}, {0,3,4,6,7,1,2,5}, {3,4,6,7,0,1,2,5},
		{0,1,2,4,6,7,3,5}, {1,2,4,6,7,0,3,5}, {0,2,4,6,7,1,3,5}, {2,4,6,7,0,1,3,5},
		{0,1,4,6,7,2,3,5}, {1,4,6,7,0,2,3,5}, {0,4,6,7,1,2,3,5}, {4,6,7,0,1,2,3,5},
		{0,1,2,3,6,7,4,5}, {1,2,3,6,7,0,4,5}, {0,2,3,6,7,1,4,5}, {2,3,6,7,0,1,4,5},
		{0,1,3,6,7,2,4,5}, {1,3,6,7,0,2,4,5}, {0,3,6,7,1,2,4,5}, {3,6,7,0,1,2,4,5},
		{0,1,2,6,7,3,4,5}, {1,2,6,7,0,3,4,5}, {0,2,6,7,1,3,4,5}, {2,6,7,0,1,3,4,5},
		{0,1,6,7,2,3,4,5}, {1,6,7,0,2,3,4,5}, {0,6,7,1,2,3,4,5}, {6,7,0,1,2,3,4,5},
		{0,1,2,3,4,5,7,6}, {1,2,3,4,5,7,0,6}, {0,2,3,4,5,7,1,6}, {2,3,4,5,7,0,1,6},
		{0,1,3,4,5,7,2,6}, {1,3,4,5,7,0,2,6}, {0,3,4,5,7,1,2,6}, {3,4,5,7,0,1,2,6},
		{0,1,2,4,5,7,3,6}, {1,2,4,5,7,0,3,6}, {0,2,4,5,7,1,3,6}, {2,4,5,7,0,1,3,6},
		{0,1,4,5,7,2,3,6}, {1,4,5,7,0,2,3,6}, {0,4,5,7,1,2,3,6}, {4,5,7,0,1,2,3,6},
		{0,1,2,3,5,7,4,6}, {1,2,3,5,7,0,4,6}, {0,2,3,5,7,1,4,6}, {2,3,5,7,0,1,4,6},
		{0,1,3,5,7,2,4,6}, {1,3,5,7,0,2,4,6}, {0,3,5,7,1,2,4,6}, {3,5,7,0,1,2,4,6},
		{0,1,2,5,7,3,4,6}, {1,2,5,7,0,3,4,6}, {0,2,5,7,1,3,4,6}, {2,5,7,0,1,3,4,6},
		{0,1,5,7,2,3,4,6}, {1,5,7,0,2,3,4,6}, {0,5,7,1,2,3,4,6}, {5,7,0,1,2,3,4,6},
		{0,1,2,3,4,7,5,6}, {1,2,3,4,7,0,5,6}, {0,2,3,4,7,1,5,6}, {2,3,4,7,0,1,5,6},
		{0,1,3,4,7,2,5,6}, {1,3,4,7,0,2,5,6}, {0,3,4,7,1,2,5,6}, {3,4,7,0,1,2,5,6},
		{0,1,2,4,7,3,5,6}, {1,2,4,7,0,3,5,6}, {0,2,4,7,1,3,5,6}, {2,4,7,0,1,3,5,6},
		{0,1,4,7,2,3,5,6}, {1,4,7,0,2,3,5,6}, {0,4,7,1,2,3,5,6}, {4,7,0,1,2,3,5,6},
		{0,1,2,3,7,4,5,6}, {1,2,3,7,0,4,5,6}, {0,2,3,7,1,4,5,6}, {2,3,7,0,1,4,5,6},
		{0,1,3,7,2,4,5,6}, {1,3,7,0,2,4,5,6}, {0,3,7,1,2,4,5,6}, {3,7,0,1,2,4,5,6},
		{0,1,2,7,3,4,5,6}, {1,2,7,0,3,4,5,6}, {0,2,7,1,3,4,5,6}, {2,7,0,1,3,4,5,6},
		{0,1,7,2,3,4,5,6}, {1,7,0,2,3,4,5,6}, {0,7,1,2,3,4,5,6}, {7,0,1,2,3,4,5,6}};

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
		vector0 = _mm_shuffle_epi8(vector0, 
			_mm_loadl_epi64((__m128i*)(table[bitmask0 & 0x7F])));
		vector1 = _mm_shuffle_epi8(vector1, 
			_mm_loadl_epi64((__m128i*)(table[(bitmask0 >> 8) & 0x7F])));

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