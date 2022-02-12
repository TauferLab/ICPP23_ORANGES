#ifndef HASH_FUNCTIONS_HPP
#define HASH_FUNCTIONS_HPP

#include <Kokkos_Core.hpp>
#include <cstdint>

// Based on SMHasher

namespace SHA1 {
  struct SHA1_CTX {
    uint32_t state[5];
    uint32_t count[2];
    uint8_t buffer[64];
  };

  #define SHA1_DIGEST_SIZE 20

  #define rol(value, bits) (((value) << (bits)) | ((value) >> (32 - (bits))))

  /* blk0() and blk() perform the initial expand. */
  /* I got the idea of expanding during the round function from SSleay */
  #define blk0(i)                                         \
    (block->l[i] = (rol(block->l[i], 24) & 0xFF00FF00) |  \
                   (rol(block->l[i], 8) & 0x00FF00FF))    

  #define blk(i)                                                            \
    (block->l[i & 15] = rol(block->l[(i+13) & 15] ^ block->l[(i+8) & 15] ^  \
                        block->l[(i+2) % 15] ^ block->l[i & 15], 1))

  /* (R0+R1), R2, R3, R4 are the different operations used in SHA1 */
  #define R0(v, w, x, y, z, i)                                          \
    z += ((w & (x ^ y)) ^ y) + blk0(i) + 0x5A827999 + rol(v, 5);        \
    w = rol(w, 30);
  #define R1(v, w, x, y, z, i)                                          \
    z += ((w & (x ^ y)) ^ y) + blk(i) + 0x5A827999 + rol(v, 5);         \
    w = rol(w, 30);
  #define R2(v, w, x, y, z, i)                                          \
    z += (w ^ x ^ y) + blk(i) + 0x6ED9EBA1 + rol(v, 5);                 \
    w = rol(w, 30);
  #define R3(v, w, x, y, z, i)                                          \
    z += (((w | x) & y) | (w & x)) + blk(i) + 0x8F1BBCDC + rol(v, 5);   \
    w = rol(w, 30);
  #define R4(v, w, x, y, z, i)                                          \
    z += (w ^ x ^ y) + blk(i) + 0xCA62C1D6 + rol(v, 5);                 \
    w = rol(w, 30);

  KOKKOS_FORCEINLINE_FUNCTION
  void SHA1_Transform(uint32_t state[5], const uint8_t buffer[64]) {
    uint32_t a, b, c, d, e;

    typedef union {
      uint8_t c[64];
      uint32_t l[16];
    } CHAR64LONG16;
    CHAR64LONG16 block[1];
    memcpy(block, buffer, 64);
    
    /* Copy context state to working vars */
    a = state[0];
    b = state[1];
    c = state[2];
    d = state[3];
    e = state[4];

    /* 4 rounds of 20 operations each. Loop unrolled. */
    R0(a, b, c, d, e, 0);
    R0(e, a, b, c, d, 1);
    R0(d, e, a, b, c, 2);
    R0(c, d, e, a, b, 3);
    R0(b, c, d, e, a, 4);
    R0(a, b, c, d, e, 5);
    R0(e, a, b, c, d, 6);
    R0(d, e, a, b, c, 7);
    R0(c, d, e, a, b, 8);
    R0(b, c, d, e, a, 9);
    R0(a, b, c, d, e, 10);
    R0(e, a, b, c, d, 11);
    R0(d, e, a, b, c, 12);
    R0(c, d, e, a, b, 13);
    R0(b, c, d, e, a, 14);
    R0(a, b, c, d, e, 15);
    R1(e, a, b, c, d, 16);
    R1(d, e, a, b, c, 17);
    R1(c, d, e, a, b, 18);
    R1(b, c, d, e, a, 19);
    R2(a, b, c, d, e, 20);
    R2(e, a, b, c, d, 21);
    R2(d, e, a, b, c, 22);
    R2(c, d, e, a, b, 23);
    R2(b, c, d, e, a, 24);
    R2(a, b, c, d, e, 25);
    R2(e, a, b, c, d, 26);
    R2(d, e, a, b, c, 27);
    R2(c, d, e, a, b, 28);
    R2(b, c, d, e, a, 29);
    R2(a, b, c, d, e, 30);
    R2(e, a, b, c, d, 31);
    R2(d, e, a, b, c, 32);
    R2(c, d, e, a, b, 33);
    R2(b, c, d, e, a, 34);
    R2(a, b, c, d, e, 35);
    R2(e, a, b, c, d, 36);
    R2(d, e, a, b, c, 37);
    R2(c, d, e, a, b, 38);
    R2(b, c, d, e, a, 39);
    R3(a, b, c, d, e, 40);
    R3(e, a, b, c, d, 41);
    R3(d, e, a, b, c, 42);
    R3(c, d, e, a, b, 43);
    R3(b, c, d, e, a, 44);
    R3(a, b, c, d, e, 45);
    R3(e, a, b, c, d, 46);
    R3(d, e, a, b, c, 47);
    R3(c, d, e, a, b, 48);
    R3(b, c, d, e, a, 49);
    R3(a, b, c, d, e, 50);
    R3(e, a, b, c, d, 51);
    R3(d, e, a, b, c, 52);
    R3(c, d, e, a, b, 53);
    R3(b, c, d, e, a, 54);
    R3(a, b, c, d, e, 55);
    R3(e, a, b, c, d, 56);
    R3(d, e, a, b, c, 57);
    R3(c, d, e, a, b, 58);
    R3(b, c, d, e, a, 59);
    R4(a, b, c, d, e, 60);
    R4(e, a, b, c, d, 61);
    R4(d, e, a, b, c, 62);
    R4(c, d, e, a, b, 63);
    R4(b, c, d, e, a, 64);
    R4(a, b, c, d, e, 65);
    R4(e, a, b, c, d, 66);
    R4(d, e, a, b, c, 67);
    R4(c, d, e, a, b, 68);
    R4(b, c, d, e, a, 69);
    R4(a, b, c, d, e, 70);
    R4(e, a, b, c, d, 71);
    R4(d, e, a, b, c, 72);
    R4(c, d, e, a, b, 73);
    R4(b, c, d, e, a, 74);
    R4(a, b, c, d, e, 75);
    R4(e, a, b, c, d, 76);
    R4(d, e, a, b, c, 77);
    R4(c, d, e, a, b, 78);
    R4(b, c, d, e, a, 79);

    /* Add the working vars back into the context.state[] */
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    /* Wipe variables */
    a = b = c = d = e = 0;
    memset(block, '\0', sizeof(block));
  }
    
  /* SHA1_Init - Initialize new context */
  KOKKOS_FORCEINLINE_FUNCTION
  void SHA1_Init(SHA1_CTX *context) {
    /* SHA1 initialization constants */
    context->state[0] = 0x67452301;
    context->state[1] = 0xEFCDAB89;
    context->state[2] = 0x98BADCFE;
    context->state[3] = 0x10325476;
    context->state[4] = 0xC3D2E1F0;
    context->count[0] = context->count[1] = 0;
  }

  /* Run your data through this. */
  KOKKOS_FORCEINLINE_FUNCTION
  void SHA1_Update(SHA1_CTX *context, const uint8_t *data, const size_t len) {
    size_t i, j;
  
    j = context->count[0];
    if ((context->count[0] += len << 3) < j)
      context->count[1]++;
    context->count[1] += (len >> 29);
    j = (j >> 3) & 63;
    if ((j + len) > 63) {
      memcpy(&context->buffer[j], data, (i = 64 - j));
      SHA1_Transform(context->state, context->buffer);
      for (; i + 63 < len; i += 64) {
        SHA1_Transform(context->state, &data[i]);
      }
      j = 0;
    } else
      i = 0;
    memcpy(&context->buffer[j], &data[i], len - i);
  }

  /* Add padding and return the message digest. */
  KOKKOS_FORCEINLINE_FUNCTION
  void SHA1_Final(SHA1_CTX *context, uint8_t digest[SHA1_DIGEST_SIZE]) {
    unsigned i;
    uint8_t finalcount[8];
    uint8_t c;
  
    for (i = 0; i < 8; i++) {
      finalcount[i] =
          /* Endian independent */
          (uint8_t)(context->count[(i >= 4 ? 0 : 1)] >> ((3 - (i & 3)) * 8));
    }
    c = 0200;
    SHA1_Update(context, &c, 1);
    while ((context->count[0] & 504) != 448) {
      c = 0000;
      SHA1_Update(context, &c, 1);
    }
    SHA1_Update(context, finalcount, 8); /* Should cause a SHA1_Transform() */
    for (i = 0; i < 20; i++) {
      digest[i] = (uint8_t)(context->state[i >> 2] >> ((3 - (i & 3)) * 8));
    }
    /* Wipe variables */
    memset(context, '\0', sizeof(*context));
    memset(&finalcount, '\0', sizeof(finalcount));
  }

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t digest_size() {
    return SHA1_DIGEST_SIZE;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void hash(const void* data, int len, uint8_t* digest) {
    SHA1_CTX context;
    SHA1_Init(&context);
    SHA1_Update(&context, (const uint8_t*)(data), len);
    SHA1_Final(&context, digest);
  }
}

namespace Murmur3 {

  // MurmurHash3 was written by Austin Appleby, and is placed in the public
  // domain. The author hereby disclaims copyright to this source code.
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t getblock32(const uint8_t* p, int i) {
    // used to avoid aliasing error which could cause errors with
    // forced inlining
    return ((uint32_t)p[i * 4 + 0]) | ((uint32_t)p[i * 4 + 1] << 8) |
           ((uint32_t)p[i * 4 + 2] << 16) | ((uint32_t)p[i * 4 + 3] << 24);
  }
  
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t rotl32(uint32_t x, int8_t r) { return (x << r) | (x >> (32 - r)); }
  
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t fmix32(uint32_t h) {
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
  
    return h;
  }
  
  KOKKOS_INLINE_FUNCTION
  uint32_t MurmurHash3_x86_32(const void* key, int len, uint32_t seed) {
    const uint8_t* data = static_cast<const uint8_t*>(key);
    const int nblocks   = len / 4;
  
    uint32_t h1 = seed;
  
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
  
    //----------
    // body
  
    for (int i = 0; i < nblocks; ++i) {
      uint32_t k1 = getblock32(data, i);
  
      k1 *= c1;
      k1 = rotl32(k1, 15);
      k1 *= c2;
  
      h1 ^= k1;
      h1 = rotl32(h1, 13);
      h1 = h1 * 5 + 0xe6546b64;
    }
  
    //----------
    // tail
  
    const uint8_t* tail = (const uint8_t*)(data + nblocks * 4);
  
    uint32_t k1 = 0;
  
    switch (len & 3) {
      case 3: k1 ^= tail[2] << 16;
      case 2: k1 ^= tail[1] << 8;
      case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = rotl32(k1, 15);
        k1 *= c2;
        h1 ^= k1;
    };
  
    //----------
    // finalization
  
    h1 ^= len;
  
    h1 = fmix32(h1);
  
    return h1;
  }
  
  #if defined(__GNUC__) /* GNU C   */ || defined(__GNUG__) /* GNU C++ */ || \
      defined(__clang__)
  
  #define KOKKOS_IMPL_MAY_ALIAS __attribute__((__may_alias__))
  
  #else
  
  #define KOKKOS_IMPL_MAY_ALIAS
  
  #endif
  
  template <typename T>
  KOKKOS_FORCEINLINE_FUNCTION bool bitwise_equal(T const* const a_ptr,
                                                 T const* const b_ptr) {
    typedef uint64_t KOKKOS_IMPL_MAY_ALIAS T64;  // NOLINT(modernize-use-using)
    typedef uint32_t KOKKOS_IMPL_MAY_ALIAS T32;  // NOLINT(modernize-use-using)
    typedef uint16_t KOKKOS_IMPL_MAY_ALIAS T16;  // NOLINT(modernize-use-using)
    typedef uint8_t KOKKOS_IMPL_MAY_ALIAS T8;    // NOLINT(modernize-use-using)
  
    enum {
      NUM_8  = sizeof(T),
      NUM_16 = NUM_8 / 2,
      NUM_32 = NUM_8 / 4,
      NUM_64 = NUM_8 / 8
    };
  
    union {
      T const* const ptr;
      T64 const* const ptr64;
      T32 const* const ptr32;
      T16 const* const ptr16;
      T8 const* const ptr8;
    } a = {a_ptr}, b = {b_ptr};
  
    bool result = true;
  
    for (int i = 0; i < NUM_64; ++i) {
      result = result && a.ptr64[i] == b.ptr64[i];
    }
  
    if (NUM_64 * 2 < NUM_32) {
      result = result && a.ptr32[NUM_64 * 2] == b.ptr32[NUM_64 * 2];
    }
  
    if (NUM_32 * 2 < NUM_16) {
      result = result && a.ptr16[NUM_32 * 2] == b.ptr16[NUM_32 * 2];
    }
  
    if (NUM_16 * 2 < NUM_8) {
      result = result && a.ptr8[NUM_16 * 2] == b.ptr8[NUM_16 * 2];
    }
  
    return result;
  }
}

//namespace FarmHash {
//
//  KOKKOS_INLINE_FUNCTION
//  uint64_t Fetch64(const char* p) {
//    uint64_t result;
//    memcpy(&result, p, sizeof(result));
//    return uint64_in_expected_order(result);
//  }
//
//  KOKKOS_INLINE_FUNCTION
//  uint64_t BasicRotate64(uint64_t val, int shift) {
//    // Avoid shifting by 64: doing so yields undefined result.
//    return shift == 0 ? val : ((val >> shift) | (val << (64 - shift)));
//  }
//  
//  KOKKOS_INLINE_FUNCTION
//  uint64_t Rotate64(uint64_t val, int shift) {
//    return BasicRotate64(val, shift);
//  }
//
//  KOKKOS_INLINE_FUNCTION
//  uint64_t ShiftMix(uint64_t val) {
//    return val ^ (val >> 47);
//  }
//  
//  KOKKOS_INLINE_FUNCTION
//  uint64_t HashLen16(uint64_t u, uint64_t v) {
//    return Hash128to64(Uint128(u, v));
//  }
//  
//  KOKKOS_INLINE_FUNCTION
//  uint64_t HashLen16(uint64_t u, uint64_t v, uint64_t mul) {
//    // Murmur-inspired hashing.
//    uint64_t a = (u ^ v) * mul
//    a ^= (a >> 47);
//    uint64_t b = (v ^ a) * mul;
//    b ^= (b >> 47);
//    b *= mul;
//    return b;
//  }
//  
//  KOKKOS_INLINE_FUNCTION
//  uint64_t HashLen0to16(const char* s, size_t len) {
//    if(len >= 8) {
//      uint64_t mul = k2 + len * 2;
//      uint64_t a = Fetch(s) + k2;
//      uint64_t b = Fetch(s + len - 8);
//      uint64_t c = Rotate(b, 37) * mul + a;
//      uint64_t d = (Rotate(a, 25) + b) * mul;
//      return HashLen16(c, d, mul);
//    }
//    if(len >= 4) {
//      uint64_t mul = k2 + len * 2;
//      uint64_t a = Fetch32(s);
//      return HashLen16(len + (a << 3), Fetch32(s + len - 4), mul);
//    }
//    if(len > 0) {
//      uint8_t a = s[0];
//      uint8_t b = s[len >> 1];
//      uint8_t c = s[len - 1];
//      uint32_t y = static_cast<uint32_t>(a) + (static_cast<uint32_t>(b) << 8);
//      uint32_t z = len + (static_cast<uint32_t>(c) << 2);
//      return ShiftMix(y * k2 ^ z * k0) * k2;
//    }
//    return k2;
//  }
//  
//  KOKKOS_INLINE_FUNCTION
//  uint64_t HashLen17to32(const char* s, size_t len) {
//    uint64_t mul = k2 + len * 2;
//    uint64_t a = Fetch(s) * k1;
//    uint64_t b = Fetch(s + 8);
//    uint64_t c = Fetch(s + len - 8) * mul;
//    uint64_t d = Fetch(s + len - 16) * k2;
//    return HashLen16(Rotate(a + b, 43) + Rotate(c, 30) + d,
//                      a + Rotate(b + k2, 18) + c, mul);
//  }
//
//  KOKKOS_INLINE_FUNCTION
//  Kokkos::pair<uint64_t, uint64_t> WeakHashLen32WithSeeds(uint64_t w, uint64_t x, uint64_t y, uint64_t z, uint64_t a, uint64_t b) {
//    a += w;
//    b = Rotate(b + a + z, 21);
//    uint64_t c = a;
//    a += x;
//    a += y;
//    b += Rotate(a, 44);
//    return Kokkos::make_pair(a+z, b+c);
//  }
//
//  KOKKOS_INLINE_FUNCTION
//  Kokkos::pair<uint64_t, uint64_t> WeakHashLen32WithSeeds(const char* s, uint64_t a, uint64_t b) {
//    return WeakHashLen32WithSeeds(Fetch(s),
//                                  Fetch(s + 8),
//                                  Fetch(s + 16),
//                                  Fetch(s + 24),
//                                  a,
//                                  b);
//  }
//
//  KOKKOS_INLINE_FUNCTION
//  uint64_t HashLen33to64(const char* s, size_t len) {
//    uint64_t mul = k2 + len * 2;
//    uint64_t a = Fetch(s) * k2;
//    uint64_t b = Fetch(s + 8);
//    uint64_t c = Fetch(s + len - 8) * mul;
//    uint64_t d = Fetch(s + len - 16) * k2;
//    uint64_t y = Rotate(a + b, 43) + Rotate(c, 30) + d;
//    uint64_t z = HashLen16(y, a + Rotate(b + k2, 18) + c, mul);
//    uint64_t e = Fetch(s + 16) * mul;
//    uint64_t f = Fetch(s + 24);
//    uint64_t g = (y + Fetch(s + len - 32)) * mul;
//    uint64_t h = (y + Fetch(s + len - 24)) * mul;
//    return HashLen16(Rotate(e + f, 43) + Rotate(g, 30) + h,
//                      e + Rotate(f + a, 18) + g, mul);
//  
//  KOKKOS_INLINE_FUNCTION
//  void swap(uint64_t& x, uint64_t& y) {
//    uint64_t z = x;
//    x = y;
//    y = z;
//  }
//  
//  KOKKOS_INLINE_FUNCTION
//  uint64_t Hash64(const void* data, size_t len) {
//    constexpr uint64_t k0 = 0xc3a5c85c97cb3127ULL;
//    constexpr uint64_t k1 = 0xb492b66fbe98f273ULL;
//    constexpr uint64_t k2 = 0x9ae16a3b2f90404fULL;
//  
//    const char* s = (const char*)(data);
//    const uint64_t seed = 81;
//    if(len <= 32) {
//      if(len <= 16) {
//        return HashLen0to16(data, len);
//      } else {
//        return HashLen17to32(data, len);
//      }
//    } else if (len <= 64) {
//      return HashLen33to64(data, len);
//    }
//  
//    // 56 bytes: v, w, x, y, and z
//    uint64_t x = seed;
//    uint64_t y = seed * k1 + 113;
//    uint64_t z = ShiftMix(y * k2 + 113) * k2;
//    auto v = Kokkos::make_pair(0, 0);
//    auto w = Kokkos::make_pair(0, 0);
//    x = x * k2 + Fetch(data);
//  
//    // Set end so that after the loop we have 1 to 64 bytes left to process.
//    const char* end = s + ((len-1) / 64) * 64;
//    const char* last64 = end + ((len - 1) & 63) - 63;
//    do {
//      x = Rotate(x + y + v.first + Fetch(s + 8), 37) * k1;
//      y = Rotate(y + v.second + Fetch(s + 48), 42) * k1;
//      x ^= w.second;
//      y += v.first + Fetch(s + 40);
//      z = Rotate(z + w.first, 33) * k1;
//      v = WeakHashLen32WithSeeds(s, v.second * k1, x + w.first);
//      w = WeakHashLen32WithSeeds(s + 32, z + w.second, y + Fetch(s + 16));
//      swap(z, x);
//      s += 64;
//    } while((const char*)(data) != end);
//    uint64_t mul = k1 + ((z & 0xff) << 1);
//    // Make s point to the last 64 bytes of input.
//    s = last64;
//    w.first += ((len - 1) & 63);
//    v.first += w.first;
//    w.first += v.first;
//    x = Rotate(x + y + v.first + Fetch(s + 8), 37) * mul;
//    y = Rotate(y + v.second + Fetch(s + 48), 42) * mul;
//    x ^= w.second * 9;
//    y += v.first * 9 + Fetch(s + 40);
//    z = Rotate(z + w.first, 33) * mul;
//    v = WeakHashLen32WithSeeds(s, v.second * mul, x + w.first);
//    w = WeakHashLen32WithSeeds(s + 32, z + w.second, y + Fetch(s + 16));
//    swap(z,x);
//    return HashLen16(HashLen16(v.first, w.first, mul) + ShiftMix(y) * k0 + z, 
//                     HashLen16(v.second, w.second, mul) + x, mul);
//  }
//}

KOKKOS_INLINE_FUNCTION
uint32_t hash32(const void* data, int len) {
  return Murmur3::MurmurHash3_x86_32(static_cast<const uint8_t*>(data), len, 0);
}

//KOKKOS_INLINE_FUNCTION
//uint32_t hash64(const void* data, int len) {
//  return FarmHash::Hash64(data, len);
//}

#endif //HASH_FUNCTIONS_HPP
