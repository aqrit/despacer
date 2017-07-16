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


size_t despace_ssse3_detect( void* dst_void, void* src_void, size_t length )
{
	uint8_t* dst = (uint8_t*)dst_void;
	uint8_t* src = (uint8_t*)src_void;
	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i mask_20 = _mm_set1_epi8( 0x20 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);

	for( ; length >= 16; length-=16 ){
		__m128i v = _mm_loadu_si128((__m128i*)src);
		int m = _mm_movemask_epi8(
			_mm_or_si128(_mm_cmpeq_epi8(mask_20, v),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v))));

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

		// p = just shift and extra bit on the last pass
		t = _mm_srli_epi64(s, 33);
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


size_t despace_ssse3_lut_512b( void* dst_void, void* src_void, size_t length )
{
	__m128i * src = (__m128i *)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	static const uint32_t table[128] __attribute__((aligned(64))) = {
		0x73625140, 0x04736251, 0x14736250, 0x15047362, 0x24736150, 0x25047361, 0x25147360, 0x26150473,
		0x34726150, 0x35047261, 0x35147260, 0x36150472, 0x35247160, 0x36250471, 0x36251470, 0x37261504,
		0x43726150, 0x45037261, 0x45137260, 0x46150372, 0x45237160, 0x46250371, 0x46251370, 0x47261503,
		0x45327160, 0x46350271, 0x46351270, 0x47361502, 0x46352170, 0x47362501, 0x47362510, 0x40372615,
		0x53726140, 0x54037261, 0x54137260, 0x56140372, 0x54237160, 0x56240371, 0x56241370, 0x57261403,
		0x54327160, 0x56340271, 0x56341270, 0x57361402, 0x56342170, 0x57362401, 0x57362410, 0x50372614,
		0x53427160, 0x56430271, 0x56431270, 0x57461302, 0x56432170, 0x57462301, 0x57462310, 0x50472613,
		0x56423170, 0x57463201, 0x57463210, 0x50473612, 0x57463120, 0x50473621, 0x51473620, 0x51403726,
		0x63725140, 0x64037251, 0x64137250, 0x65140372, 0x64237150, 0x65240371, 0x65241370, 0x67251403,
		0x64327150, 0x65340271, 0x65341270, 0x67351402, 0x65342170, 0x67352401, 0x67352410, 0x60372514,
		0x63427150, 0x65430271, 0x65431270, 0x67451302, 0x65432170, 0x67452301, 0x67452310, 0x60472513,
		0x65423170, 0x67453201, 0x67453210, 0x60473512, 0x67453120, 0x60473521, 0x61473520, 0x61403725,
		0x63527140, 0x64530271, 0x64531270, 0x67541302, 0x64532170, 0x67542301, 0x67542310, 0x60572413,
		0x64523170, 0x67543201, 0x67543210, 0x60573412, 0x67543120, 0x60573421, 0x61573420, 0x61503724,
		0x63524170, 0x67534201, 0x67534210, 0x60574312, 0x67534120, 0x60574321, 0x61574320, 0x61504723,
		0x67524130, 0x60574231, 0x61574230, 0x61504732, 0x62574130, 0x62504731, 0x62514730, 0x62514037};

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

		// sort
		int bitmask0 = _mm_movemask_epi8(bytemask0);
		uint64_t shuf0 = table[bitmask0 & 0x7F];
		uint64_t shuf1 = table[(bitmask0 >> 8) & 0x7F];
		shuf0 |= shuf0 << 28; // unpack
		shuf1 |= shuf1 << 28;
		vector0 = _mm_shuffle_epi8( vector0, _mm_cvtsi64_si128( shuf0 ) );
		vector1 = _mm_shuffle_epi8( vector1, _mm_cvtsi64_si128( shuf1 ) );

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

// embed the popcnt directly into the shuffle table
//
// slower than _mm_sad_epu8 (here) 
// because waiting on load ? and/or
// pshufb xmm,m128 is faster than load + pshufb xmm,xmm ?
__m128i lut_1mb_2 [0x10000];
void gen_table_1mb_2( void ){
	for( int i = 0; i < 0x10000; i++ ){
		uint8_t* p = (uint8_t*) &lut_1mb_2[i];
		p[15] = 16 - __builtin_popcount(i);
		for( size_t j = 0; j < 16; j++ ){
			if( !(i & (1 << j)) ){
				*p++ = j+1;
			}
		}
	}
}
size_t despace_ssse3_lut_1mb_2( void* dst_void, void* src_void, size_t length )
{
	__m128i * src = (__m128i *)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i mask_20 = _mm_set1_epi8( 0x20 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);
	const __m128i mask_FF = _mm_cmpeq_epi8(_mm_setzero_si128(), _mm_setzero_si128());

	for( ; length >= 16; length-=16 ){
		__m128i v = _mm_loadu_si128(src++);
		__m128i bytemask = _mm_or_si128(_mm_cmpeq_epi8(mask_20, v),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v)));
		bytemask = _mm_load_si128(&lut_1mb_2[_mm_movemask_epi8(bytemask)]);
		_mm_storeu_si128((__m128i*)dst, _mm_shuffle_epi8(v, _mm_add_epi8(bytemask, mask_FF)));
		dst += _mm_cvtsi128_si32(_mm_bsrli_si128(bytemask, 15));
	}
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}
