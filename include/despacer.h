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
	__m128i * src = (__m128i *)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m128i mask_70 = _mm_set1_epi8( 0x70 );
	const __m128i mask_20 = _mm_set1_epi8( 0x20 );
	const __m128i lut_cntrl = _mm_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0xFF, 0xFF, 0x00, 0x00, 0xFF, 0x00, 0x00);

	for( ; length >= 16; length-=16 ){
		__m128i v = _mm_loadu_si128(src++);
		__m128i bytemask = _mm_or_si128(_mm_cmpeq_epi8(mask_20, v),
			_mm_shuffle_epi8(lut_cntrl, _mm_adds_epu8(mask_70, v)));
		int bitmask = _mm_movemask_epi8(bytemask);
		v = _mm_shuffle_epi8(v, _mm_load_si128(&lut_1mb[bitmask]));
		_mm_storeu_si128((__m128i*)dst, v);
		__m128i hsum = _mm_sad_epu8(_mm_setzero_si128(), bytemask);
		hsum = _mm_add_epi32(hsum, _mm_unpackhi_epi64(hsum, hsum));
		dst += (uint8_t)(16 + _mm_cvtsi128_si32(hsum));
	}
	dst += despace_branchless(dst, src, length);
	return (size_t)(dst - ((uint8_t*)dst_void));
}
