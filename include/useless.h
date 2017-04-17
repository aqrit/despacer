// despacers with mediocre results [for posterity's sake]
//
//


// force use of setcc instruction (slow)
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


// shift block instead of re-reading byte
size_t despace_block_shift( void* dst_void, void* src_void, size_t length )
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
			*dst = (uint8_t)asrc;
			asrc >>= 8;
			dst += mask & 1;
			mask >>= 8;
		}
		src += 8;
	}
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


// pcmpestrm seems slower than the sse3 detection but uses fewer registers
size_t despace_sse42_detect( void* dst_void, void* src_void, size_t length )
{
	uint8_t* dst = (uint8_t*)dst_void;
	uint8_t* src = (uint8_t*)src_void;
	const __m128i targetchars = _mm_cvtsi32_si128( 0x090A0D20 );

	for( ; length >= 16; length-=16 ){
		int m = _mm_cvtsi128_si32( _mm_cmpestrm( targetchars, 4, _mm_loadu_si128((__m128i*)src), 16, 
				_SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_ANY | _SIDD_BIT_MASK | _SIDD_NEGATIVE_POLARITY) );

		for( int i = 0; i < 16; i++ ){
			*dst = src[i];
			dst += m & 1;
			m >>= 1;
		}
		src += 16;
	}
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


// could write full block instead of memmov if dst and src didn't overlap...
// _mm_cmpistri might be faster?
size_t despace_sse42_scan( void* dst_void, void* src_void, size_t length )
{
	uint8_t* dst = (uint8_t*)dst_void;
	uint8_t* src = (uint8_t*)src_void;
	const __m128i targetchars = _mm_cvtsi32_si128( 0x090A0D20 );

	size_t pos = 0;
	while( pos <= length-16 ){
		int r = _mm_cmpestri( targetchars, 4, _mm_loadu_si128((__m128i*)&src[pos]), 16, 0 );
		memmove(dst, &src[pos], r);
		pos += r;
		dst += r;
		if( r != 16 ) ++pos; // skip space if any
	}
	dst += despace_branchless(dst, &src[pos], length - pos);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


// slower than just creating mask with sse2 then using scalar
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
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}


size_t despace_ssse3_cumsum( void* dst_void, void* src_void, size_t length )
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
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}
