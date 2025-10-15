/**
 * SNOVA 2025: Combination of OQS AES functions in a single file
 *
 * Edits: Declared internal functions as static.
 *
 * Copyright: see below
 */

#include "snova.h"
#define AES128_CTR SNOVA_NAMESPACE(AES128_CTR)
#define AES256_ECB SNOVA_NAMESPACE(AES256_ECB)

#ifdef OPENSSL_AES

// NIST rng.c

#include <openssl/err.h>
#include <openssl/evp.h>
#include <string.h>

void handleErrors(void) {
	ERR_print_errors_fp(stderr);
	abort();
}

// Use whatever AES implementation you have. This uses AES from openSSL library
//    key - 256-bit AES key
//    ctr - a 128-bit plaintext value
//    buffer - a 128-bit ciphertext value
void AES256_ECB(unsigned char *key, unsigned char *ctr, unsigned char *buffer) {
	EVP_CIPHER_CTX *ctx;

	int len;

	/* Create and initialise the context */
	if (!(ctx = EVP_CIPHER_CTX_new())) {
		handleErrors();
	}

	if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_ecb(), NULL, key, NULL)) {
		handleErrors();
	}

	if (1 != EVP_EncryptUpdate(ctx, buffer, &len, ctr, 16)) {
		handleErrors();
	}

	/* Clean up */
	EVP_CIPHER_CTX_free(ctx);
}

int AES128_CTR(unsigned char *data, size_t num_bytes, const unsigned char *input, size_t inputByteLen) {
	(void)inputByteLen;
	EVP_CIPHER_CTX *context;
	int len;

	memset(data, 0, num_bytes);
	context = EVP_CIPHER_CTX_new();
	EVP_EncryptInit_ex(context, EVP_aes_128_ctr(), NULL, input, NULL);
	EVP_EncryptUpdate(context, data, &len, data, num_bytes);
	EVP_CIPHER_CTX_free(context);

	return 0;
}

#else

#pragma GCC diagnostic ignored "-Wunused-function"

/**
 * \file aes.h
 * \brief Header defining the API for OQS AES; not part of the public OQS API
 *
 * <b>Note this is not part of the OQS public API: implementations within liboqs can use these
 * functions, but external consumers of liboqs should not use these functions.</b>
 *
 * SPDX-License-Identifier: MIT
 */

#ifndef OQS_AES_H
#define OQS_AES_H

#include <stdint.h>
#include <stdlib.h>

// #include "oqs_common.h"
/**
 * \file common.h
 * \brief Utility functions for use in liboqs.
 *
 * SPDX-License-Identifier: MIT
 */

#ifndef OQS_COMMON_H
#define OQS_COMMON_H

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Macro for terminating the program if x is
 * a null pointer.
 */
#define OQS_EXIT_IF_NULLPTR(x, loc)                                                   \
    do {                                                                              \
        if ((x) == (void *)0) {                                                       \
            fprintf(stderr, "Unexpected NULL returned from %s API. Exiting.\n", loc); \
            exit(EXIT_FAILURE);                                                       \
        }                                                                             \
    } while (0)

/**
 * This macro is intended to replace those assert()s
 * involving side-effecting statements in aes/aes_ossl.c.
 *
 * assert() becomes a no-op when -DNDEBUG is defined,
 * which causes compilation failures when the statement
 * being checked also results in side-effects.
 *
 * This is a temporary workaround until a better error
 * handling strategy is developed.
 */
#define OQS_OPENSSL_GUARD(x)                                                           \
    do {                                                                               \
        if (1 != (x)) {                                                                \
            fprintf(stderr, "Error return value from OpenSSL API: %d. Exiting.\n", x); \
            exit(EXIT_FAILURE);                                                        \
        }                                                                              \
    } while (0)

/**
 * Certain functions (such as OQS_randombytes_openssl in
 * src/rand/rand.c) take in a size_t parameter, but can
 * only handle values up to INT_MAX for those parameters.
 * This macro is a temporary workaround for such functions.
 */
#define SIZE_T_TO_INT_OR_EXIT(size_t_var_name, int_var_name) \
    int int_var_name = 0;                                    \
    if (size_t_var_name <= INT_MAX) {                        \
        int_var_name = (int)size_t_var_name;                 \
    } else {                                                 \
        exit(EXIT_FAILURE);                                  \
    }

/**
 * Defines which functions should be exposed outside the LibOQS library
 *
 * By default the visibility of all the symbols is defined to "hidden"
 * Only the library API should be marked as default
 *
 * Example: OQS_API return_value function_name(void);
 */
#if defined(_WIN32)
#define OQS_API __declspec(dllexport)
#else
#define OQS_API __attribute__((visibility("default")))
#endif

#if defined(OQS_SYS_UEFI)
#undef OQS_API
#define OQS_API
#endif

/**
 * Represents return values from functions.
 *
 * Callers should compare with the symbol rather than the individual value.
 * For example,
 *
 *     ret = OQS_KEM_encaps(...);
 *     if (ret == OQS_SUCCESS) { ... }
 *
 * rather than
 *
 *     if (!OQS_KEM_encaps(...) { ... }
 *
 */
typedef enum {
	/** Used to indicate that some undefined error occurred. */
	OQS_ERROR = -1,
	/** Used to indicate successful return from function. */
	OQS_SUCCESS = 0,
	/** Used to indicate failures in external libraries (e.g., OpenSSL). */
	OQS_EXTERNAL_LIB_ERROR_OPENSSL = 50,
} OQS_STATUS;

/**
 * CPU runtime detection flags
 */
typedef enum {
	OQS_CPU_EXT_INIT, /* Must be first */
	/* Start extension list */
	OQS_CPU_EXT_ADX,
	OQS_CPU_EXT_AES,
	OQS_CPU_EXT_AVX,
	OQS_CPU_EXT_AVX2,
	OQS_CPU_EXT_AVX512,
	OQS_CPU_EXT_BMI1,
	OQS_CPU_EXT_BMI2,
	OQS_CPU_EXT_PCLMULQDQ,
	OQS_CPU_EXT_VPCLMULQDQ,
	OQS_CPU_EXT_POPCNT,
	OQS_CPU_EXT_SSE,
	OQS_CPU_EXT_SSE2,
	OQS_CPU_EXT_SSE3,
	OQS_CPU_EXT_ARM_AES,
	OQS_CPU_EXT_ARM_SHA2,
	OQS_CPU_EXT_ARM_SHA3,
	OQS_CPU_EXT_ARM_NEON,
	/* End extension list */
	OQS_CPU_EXT_COUNT, /* Must be last */
} OQS_CPU_EXT;

/**
 * Zeros out `len` bytes of memory starting at `ptr`, then frees `ptr`.
 *
 * Can be called with `ptr = NULL`, in which case no operation is performed.
 *
 * Designed to be protected against optimizing compilers which try to remove
 * "unnecessary" operations.  Should be used for all buffers containing secret
 * data.
 *
 * @param[in] ptr The start of the memory to zero out and free.
 * @param[in] len The number of bytes to zero out.
 */
static void OQS_MEM_secure_free(void *ptr, size_t len) {
	if (ptr != NULL) {
		memset(ptr, 0, len);
		free(ptr);  // IGNORE free-check
	}
}

/**
 * Internal implementation of C11 aligned_alloc to work around compiler quirks.
 *
 * Allocates size bytes of uninitialized memory with a base pointer that is
 * a multiple of alignment. Alignment must be a power of two and a multiple
 * of sizeof(void *). Size must be a multiple of alignment.
 */
void *OQS_MEM_aligned_alloc(size_t alignment, size_t size);

/**
 * Free memory allocated with OQS_MEM_aligned_alloc.
 */
void OQS_MEM_aligned_free(void *ptr);

#if defined(__cplusplus)
}  // extern "C"
#endif

#endif  // OQS_COMMON_H

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Function to fill a key schedule given an initial key for use in ECB mode.
 *
 * @param key            Initial Key.
 * @param ctx            Abstract data structure for a key schedule.
 */
static void oqs_aes128_ECB_load_schedule(const uint8_t *key, void **ctx);

/**
 * Function to initialize a context and fill a key schedule given an initial key for
 * use in CTR mode.
 *
 * @param key            Initial Key.
 * @param ctx            Abstract data structure for a key schedule.
 */
static void oqs_aes128_CTR_inc_init(const uint8_t *key, void **ctx);

/**
 * Function to fill a context given an IV for use in CTR mode.
 *
 * Handles a 12- or 16-byte IV.  If a 12-byte IV is given, then 4 counter
 * bytes are initialized to all zeros.
 *
 * @param iv             Initialization Vector.
 * @param iv_len         Length of the initialization vector.
 * @param ctx            Abstract data structure for IV.
 */
static void oqs_aes128_CTR_inc_iv(const uint8_t *iv, size_t iv_len, void *ctx);

/**
 * Function to fill a context given an IV for use in CTR mode.
 * Handles an 8-byte IV passed as a 64-bit unsigned integer,
 * counter bytes are initialized to zero.
 *
 * @param iv             Initialization Vector as 64-bit integer.
 * @param ctx            Abstract data structure for IV.
 */
static void oqs_aes128_CTR_inc_ivu64(uint64_t iv, void *ctx);

/**
 * Function to free a key schedule.
 *
 * @param ctx            Context generated with OQS_AES128_ECB_load_schedule().
 */
static void oqs_aes128_free_schedule(void *ctx);

/**
 * Function to encrypt blocks of plaintext using ECB mode.
 * A schedule based on the key is generated and used internally.
 *
 * @param plaintext     Plaintext to be encrypted.
 * @param plaintext_len Length on the plaintext in bytes. Must be a multiple of 16.
 * @param key           Key to be used for encryption.
 * @param ciphertext    Pointer to a block of memory which >= in size to the plaintext block. The result will be written here.
 * @warning plaintext_len must be a multiple of 16.
 */
static void oqs_aes128_ECB_enc(const uint8_t *plaintext, const size_t plaintext_len, const uint8_t *key, uint8_t *ciphertext);

/**
 * Same as OQS_AES128_ECB_enc() except a schedule generated by
 * OQS_AES128_ECB_load_schedule() is passed rather then a key. This is faster
 * if the same schedule is used for multiple encryptions since it does
 * not have to be regenerated from the key.
 */
static void oqs_aes128_ECB_enc_sch(const uint8_t *plaintext, const size_t plaintext_len, const void *schedule,
                                   uint8_t *ciphertext);

/**
 * AES counter mode keystream generator.  A context generated by
 * OQS_AES128_CTR_inc_init() is passed rather then a key.
 *
 * Handles a 12- or 16-byte IV.  If a 12-byte IV is given, then 4 counter
 * bytes are initialized to all zeros.
 *
 * @param iv       12- or 16-byte initialization vector.
 * @param iv_len   Lengh of IV in bytes.
 * @param ctx      Abstract data structure for a key schedule.
 * @param out      Pointer to a block of memory which is big enough to contain out_len bytes; the result will be written here.
 * @param out_len  Length of output bytes to generate.
 */
static void oqs_aes128_CTR_inc_stream_iv(const uint8_t *iv, size_t iv_len, const void *ctx, uint8_t *out, size_t out_len);

/**
 * Function to fill a key schedule given an initial key for use in ECB mode encryption.
 *
 * @param key            Initial Key.
 * @param ctx            Abstract data structure for a key schedule.
 */
static void oqs_aes256_ECB_load_schedule(const uint8_t *key, void **ctx);

/**
 * Function to initialize a context and fill a key schedule given an initial key for
 * use in CTR mode.
 *
 * @param key            Initial Key.
 * @param ctx            Abstract data structure for a key schedule.
 */
static void oqs_aes256_CTR_inc_init(const uint8_t *key, void **ctx);

/**
 * Function to fill a context given an IV for use in CTR mode.
 *
 * Handles a 12- or 16-byte IV.  If a 12-byte IV is given, then 4 counter
 * bytes are initialized to all zeros.
 *
 * @param iv             Initialization Vector.
 * @param iv_len         Length of the initialization vector.
 * @param ctx            Abstract data structure for IV.
 */
static void oqs_aes256_CTR_inc_iv(const uint8_t *iv, size_t iv_len, void *ctx);

/**
 * Function to fill a context given an IV for use in CTR mode.
 * Handles an 8-byte IV passed as a 64-bit unsigned integer,
 * counter bytes are initialized to zero.
 *
 * @param iv             Initialization Vector as 64-bit integer.
 * @param ctx            Abstract data structure for IV.
 */
static void oqs_aes256_CTR_inc_ivu64(uint64_t iv, void *ctx);

/**
 * Function to free a key schedule.
 *
 * @param ctx            Schedule generated with OQS_AES256_ECB_load_schedule
 *                       or OQS_AES256_CTR_inc_init.
 */
static void oqs_aes256_free_schedule(void *ctx);

/**
 * Function to encrypt blocks of plaintext using ECB mode.
 * A schedule based on the key is generated and used internally.
 *
 * @param plaintext     Plaintext to be encrypted.
 * @param plaintext_len Length on the plaintext in bytes. Must be a multiple of 16.
 * @param key           Key to be used for encryption.
 * @param ciphertext    Pointer to a block of memory which >= in size to the plaintext block. The result will be written here.
 * @warning plaintext_len must be a multiple of 16.
 */
static void oqs_aes256_ECB_enc(const uint8_t *plaintext, const size_t plaintext_len, const uint8_t *key, uint8_t *ciphertext);

/**
 * Same as OQS_AES256_ECB_enc() except a schedule generated by
 * OQS_AES256_ECB_load_schedule() is passed rather then a key. This is faster
 * if the same schedule is used for multiple encryptions since it does
 * not have to be regenerated from the key.
 */
static void oqs_aes256_ECB_enc_sch(const uint8_t *plaintext, const size_t plaintext_len, const void *schedule,
                                   uint8_t *ciphertext);

/**
 * AES counter mode keystream generator.  A context generated by
 * OQS_AES256_CTR_inc_init() is passed rather then a key.
 *
 * Handles a 12- or 16-byte IV.  If a 12-byte IV is given, then 4 counter
 * bytes are initialized to all zeros.
 *
 * @param iv       12- or 16-byte initialization vector.
 * @param iv_len   Lengh of IV in bytes.
 * @param ctx      Abstract data structure for a key schedule.
 * @param out      Pointer to a block of memory which is big enough to contain out_len bytes; the result will be written here.
 * @param out_len  Length of output bytes to generate.
 */
static void oqs_aes256_CTR_inc_stream_iv(const uint8_t *iv, size_t iv_len, const void *ctx, uint8_t *out, size_t out_len);

/**
 * AES counter mode keystream generator. A context generated by
 * OQS_AES256_CTR_inc_init() and OQS_AES256_CTR_inc_iv() is passed
 * rather than a key and an IV. The counter is internally updated, which allows
 * the function to be called multiple times.
 *
 * @param ctx         Abstract data structure for key schedule and IV.
 * @param out         Pointer to a block of memory which is big enough to contain out_blks*16 bytes; the result will be written
 * here.
 * @param out_blks    Length of output blocks to generate, where one block is 16 bytes.
 */
static void oqs_aes256_CTR_inc_stream_blks(void *ctx, uint8_t *out, size_t out_blks);

#if defined(__cplusplus)
}  // extern "C"
#endif

#endif  // OQS_AES_H
// SPDX-License-Identifier: MIT
/* Adapted for OQS from PQClean. */
/*
 * AES implementation based on code from BearSSL (https://bearssl.org/)
 * by Thomas Pornin.
 *
 *
 * Copyright (c) 2016 Thomas Pornin <pornin@bolet.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// #include "aes.h"
// #include "oqs_common.h"

#define AES128_KEYBYTES 16
#define AES256_KEYBYTES 32
#define AESCTR_NONCEBYTES 12
#define AES_BLOCKBYTES 16

#define PQC_AES128_STATESIZE 88
typedef struct {
	uint64_t sk_exp[PQC_AES128_STATESIZE];
	uint8_t iv[AES_BLOCKBYTES];
} aes128ctx;

#define PQC_AES256_STATESIZE 120
typedef struct {
	uint64_t sk_exp[PQC_AES256_STATESIZE];
	uint8_t iv[AES_BLOCKBYTES];
} aes256ctx;

typedef struct {
	uint32_t sk_exp[60];
	uint8_t iv[16];
} aes256ctx_nobitslice;

static inline uint32_t br_dec32le(const unsigned char *src) {
	return (uint32_t)src[0] | ((uint32_t)src[1] << 8) | ((uint32_t)src[2] << 16) | ((uint32_t)src[3] << 24);
}

static void br_range_dec32le(uint32_t *v, size_t num, const unsigned char *src) {
	while (num-- > 0) {
		*v++ = br_dec32le(src);
		src += 4;
	}
}

static inline uint32_t br_swap32(uint32_t x) {
	x = ((x & (uint32_t)0x00FF00FF) << 8) | ((x >> 8) & (uint32_t)0x00FF00FF);
	return (x << 16) | (x >> 16);
}

static inline void br_enc32le(unsigned char *dst, uint32_t x) {
	dst[0] = (unsigned char)x;
	dst[1] = (unsigned char)(x >> 8);
	dst[2] = (unsigned char)(x >> 16);
	dst[3] = (unsigned char)(x >> 24);
}

static inline void br_enc32be(unsigned char *dst, uint32_t x) {
	dst[0] = (unsigned char)(x >> 24);
	dst[1] = (unsigned char)(x >> 16);
	dst[2] = (unsigned char)(x >> 8);
	dst[3] = (unsigned char)x;
}

static void br_range_enc32le(unsigned char *dst, const uint32_t *v, size_t num) {
	while (num-- > 0) {
		br_enc32le(dst, *v++);
		dst += 4;
	}
}

static void br_aes_ct64_bitslice_Sbox(uint64_t *q) {
	/*
	 * This S-box implementation is a straightforward translation of
	 * the circuit described by Boyar and Peralta in "A new
	 * combinational logic minimization technique with applications
	 * to cryptology" (https://eprint.iacr.org/2009/191.pdf).
	 *
	 * Note that variables x* (input) and s* (output) are numbered
	 * in "reverse" order (x0 is the high bit, x7 is the low bit).
	 */

	uint64_t x0, x1, x2, x3, x4, x5, x6, x7;
	uint64_t y1, y2, y3, y4, y5, y6, y7, y8, y9;
	uint64_t y10, y11, y12, y13, y14, y15, y16, y17, y18, y19;
	uint64_t y20, y21;
	uint64_t z0, z1, z2, z3, z4, z5, z6, z7, z8, z9;
	uint64_t z10, z11, z12, z13, z14, z15, z16, z17;
	uint64_t t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
	uint64_t t10, t11, t12, t13, t14, t15, t16, t17, t18, t19;
	uint64_t t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
	uint64_t t30, t31, t32, t33, t34, t35, t36, t37, t38, t39;
	uint64_t t40, t41, t42, t43, t44, t45, t46, t47, t48, t49;
	uint64_t t50, t51, t52, t53, t54, t55, t56, t57, t58, t59;
	uint64_t t60, t61, t62, t63, t64, t65, t66, t67;
	uint64_t s0, s1, s2, s3, s4, s5, s6, s7;

	x0 = q[7];
	x1 = q[6];
	x2 = q[5];
	x3 = q[4];
	x4 = q[3];
	x5 = q[2];
	x6 = q[1];
	x7 = q[0];

	/*
	 * Top linear transformation.
	 */
	y14 = x3 ^ x5;
	y13 = x0 ^ x6;
	y9 = x0 ^ x3;
	y8 = x0 ^ x5;
	t0 = x1 ^ x2;
	y1 = t0 ^ x7;
	y4 = y1 ^ x3;
	y12 = y13 ^ y14;
	y2 = y1 ^ x0;
	y5 = y1 ^ x6;
	y3 = y5 ^ y8;
	t1 = x4 ^ y12;
	y15 = t1 ^ x5;
	y20 = t1 ^ x1;
	y6 = y15 ^ x7;
	y10 = y15 ^ t0;
	y11 = y20 ^ y9;
	y7 = x7 ^ y11;
	y17 = y10 ^ y11;
	y19 = y10 ^ y8;
	y16 = t0 ^ y11;
	y21 = y13 ^ y16;
	y18 = x0 ^ y16;

	/*
	 * Non-linear section.
	 */
	t2 = y12 & y15;
	t3 = y3 & y6;
	t4 = t3 ^ t2;
	t5 = y4 & x7;
	t6 = t5 ^ t2;
	t7 = y13 & y16;
	t8 = y5 & y1;
	t9 = t8 ^ t7;
	t10 = y2 & y7;
	t11 = t10 ^ t7;
	t12 = y9 & y11;
	t13 = y14 & y17;
	t14 = t13 ^ t12;
	t15 = y8 & y10;
	t16 = t15 ^ t12;
	t17 = t4 ^ t14;
	t18 = t6 ^ t16;
	t19 = t9 ^ t14;
	t20 = t11 ^ t16;
	t21 = t17 ^ y20;
	t22 = t18 ^ y19;
	t23 = t19 ^ y21;
	t24 = t20 ^ y18;

	t25 = t21 ^ t22;
	t26 = t21 & t23;
	t27 = t24 ^ t26;
	t28 = t25 & t27;
	t29 = t28 ^ t22;
	t30 = t23 ^ t24;
	t31 = t22 ^ t26;
	t32 = t31 & t30;
	t33 = t32 ^ t24;
	t34 = t23 ^ t33;
	t35 = t27 ^ t33;
	t36 = t24 & t35;
	t37 = t36 ^ t34;
	t38 = t27 ^ t36;
	t39 = t29 & t38;
	t40 = t25 ^ t39;

	t41 = t40 ^ t37;
	t42 = t29 ^ t33;
	t43 = t29 ^ t40;
	t44 = t33 ^ t37;
	t45 = t42 ^ t41;
	z0 = t44 & y15;
	z1 = t37 & y6;
	z2 = t33 & x7;
	z3 = t43 & y16;
	z4 = t40 & y1;
	z5 = t29 & y7;
	z6 = t42 & y11;
	z7 = t45 & y17;
	z8 = t41 & y10;
	z9 = t44 & y12;
	z10 = t37 & y3;
	z11 = t33 & y4;
	z12 = t43 & y13;
	z13 = t40 & y5;
	z14 = t29 & y2;
	z15 = t42 & y9;
	z16 = t45 & y14;
	z17 = t41 & y8;

	/*
	 * Bottom linear transformation.
	 */
	t46 = z15 ^ z16;
	t47 = z10 ^ z11;
	t48 = z5 ^ z13;
	t49 = z9 ^ z10;
	t50 = z2 ^ z12;
	t51 = z2 ^ z5;
	t52 = z7 ^ z8;
	t53 = z0 ^ z3;
	t54 = z6 ^ z7;
	t55 = z16 ^ z17;
	t56 = z12 ^ t48;
	t57 = t50 ^ t53;
	t58 = z4 ^ t46;
	t59 = z3 ^ t54;
	t60 = t46 ^ t57;
	t61 = z14 ^ t57;
	t62 = t52 ^ t58;
	t63 = t49 ^ t58;
	t64 = z4 ^ t59;
	t65 = t61 ^ t62;
	t66 = z1 ^ t63;
	s0 = t59 ^ t63;
	s6 = t56 ^ ~t62;
	s7 = t48 ^ ~t60;
	t67 = t64 ^ t65;
	s3 = t53 ^ t66;
	s4 = t51 ^ t66;
	s5 = t47 ^ t65;
	s1 = t64 ^ ~s3;
	s2 = t55 ^ ~t67;

	q[7] = s0;
	q[6] = s1;
	q[5] = s2;
	q[4] = s3;
	q[3] = s4;
	q[2] = s5;
	q[1] = s6;
	q[0] = s7;
}

static void br_aes_ct64_ortho(uint64_t *q) {
#define SWAPN(cl, ch, s, x, y)                                      \
    do {                                                            \
        uint64_t a, b;                                              \
        a = (x);                                                    \
        b = (y);                                                    \
        (x) = (a & (uint64_t)(cl)) | ((b & (uint64_t)(cl)) << (s)); \
        (y) = ((a & (uint64_t)(ch)) >> (s)) | (b & (uint64_t)(ch)); \
    } while (0)

#define SWAP2(x, y) SWAPN(0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 1, x, y)
#define SWAP4(x, y) SWAPN(0x3333333333333333, 0xCCCCCCCCCCCCCCCC, 2, x, y)
#define SWAP8(x, y) SWAPN(0x0F0F0F0F0F0F0F0F, 0xF0F0F0F0F0F0F0F0, 4, x, y)

	SWAP2(q[0], q[1]);
	SWAP2(q[2], q[3]);
	SWAP2(q[4], q[5]);
	SWAP2(q[6], q[7]);

	SWAP4(q[0], q[2]);
	SWAP4(q[1], q[3]);
	SWAP4(q[4], q[6]);
	SWAP4(q[5], q[7]);

	SWAP8(q[0], q[4]);
	SWAP8(q[1], q[5]);
	SWAP8(q[2], q[6]);
	SWAP8(q[3], q[7]);
}

static void br_aes_ct64_interleave_in(uint64_t *q0, uint64_t *q1, const uint32_t *w) {
	uint64_t x0, x1, x2, x3;

	x0 = w[0];
	x1 = w[1];
	x2 = w[2];
	x3 = w[3];
	x0 |= (x0 << 16);
	x1 |= (x1 << 16);
	x2 |= (x2 << 16);
	x3 |= (x3 << 16);
	x0 &= (uint64_t)0x0000FFFF0000FFFF;
	x1 &= (uint64_t)0x0000FFFF0000FFFF;
	x2 &= (uint64_t)0x0000FFFF0000FFFF;
	x3 &= (uint64_t)0x0000FFFF0000FFFF;
	x0 |= (x0 << 8);
	x1 |= (x1 << 8);
	x2 |= (x2 << 8);
	x3 |= (x3 << 8);
	x0 &= (uint64_t)0x00FF00FF00FF00FF;
	x1 &= (uint64_t)0x00FF00FF00FF00FF;
	x2 &= (uint64_t)0x00FF00FF00FF00FF;
	x3 &= (uint64_t)0x00FF00FF00FF00FF;
	*q0 = x0 | (x2 << 8);
	*q1 = x1 | (x3 << 8);
}

static void br_aes_ct64_interleave_out(uint32_t *w, uint64_t q0, uint64_t q1) {
	uint64_t x0, x1, x2, x3;

	x0 = q0 & (uint64_t)0x00FF00FF00FF00FF;
	x1 = q1 & (uint64_t)0x00FF00FF00FF00FF;
	x2 = (q0 >> 8) & (uint64_t)0x00FF00FF00FF00FF;
	x3 = (q1 >> 8) & (uint64_t)0x00FF00FF00FF00FF;
	x0 |= (x0 >> 8);
	x1 |= (x1 >> 8);
	x2 |= (x2 >> 8);
	x3 |= (x3 >> 8);
	x0 &= (uint64_t)0x0000FFFF0000FFFF;
	x1 &= (uint64_t)0x0000FFFF0000FFFF;
	x2 &= (uint64_t)0x0000FFFF0000FFFF;
	x3 &= (uint64_t)0x0000FFFF0000FFFF;
	w[0] = (uint32_t)x0 | (uint32_t)(x0 >> 16);
	w[1] = (uint32_t)x1 | (uint32_t)(x1 >> 16);
	w[2] = (uint32_t)x2 | (uint32_t)(x2 >> 16);
	w[3] = (uint32_t)x3 | (uint32_t)(x3 >> 16);
}

static const unsigned char Rcon[] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36};

static uint32_t sub_word(uint32_t x) {
	uint64_t q[8];

	memset(q, 0, sizeof q);
	q[0] = x;
	br_aes_ct64_ortho(q);
	br_aes_ct64_bitslice_Sbox(q);
	br_aes_ct64_ortho(q);
	return (uint32_t)q[0];
}

static void br_aes_ct64_keysched(uint64_t *comp_skey, const unsigned char *key, unsigned int key_len) {
	unsigned int i, j, k, nk, nkf;
	uint32_t tmp;
	uint32_t skey[60];
	unsigned nrounds = 10 + ((key_len - 16) >> 2);

	nk = (key_len >> 2);
	nkf = ((nrounds + 1) << 2);
	br_range_dec32le(skey, (key_len >> 2), key);
	tmp = skey[(key_len >> 2) - 1];
	for (i = nk, j = 0, k = 0; i < nkf; i++) {
		if (j == 0) {
			tmp = (tmp << 24) | (tmp >> 8);
			tmp = sub_word(tmp) ^ Rcon[k];
		} else if (nk > 6 && j == 4) {
			tmp = sub_word(tmp);
		}
		tmp ^= skey[i - nk];
		skey[i] = tmp;
		if (++j == nk) {
			j = 0;
			k++;
		}
	}

	for (i = 0, j = 0; i < nkf; i += 4, j += 2) {
		uint64_t q[8];

		br_aes_ct64_interleave_in(&q[0], &q[4], skey + i);
		q[1] = q[0];
		q[2] = q[0];
		q[3] = q[0];
		q[5] = q[4];
		q[6] = q[4];
		q[7] = q[4];
		br_aes_ct64_ortho(q);
		comp_skey[j + 0] = (q[0] & (uint64_t)0x1111111111111111) | (q[1] & (uint64_t)0x2222222222222222) |
		                   (q[2] & (uint64_t)0x4444444444444444) | (q[3] & (uint64_t)0x8888888888888888);
		comp_skey[j + 1] = (q[4] & (uint64_t)0x1111111111111111) | (q[5] & (uint64_t)0x2222222222222222) |
		                   (q[6] & (uint64_t)0x4444444444444444) | (q[7] & (uint64_t)0x8888888888888888);
	}
}

static void br_aes_ct64_skey_expand(uint64_t *skey, const uint64_t *comp_skey, unsigned int nrounds) {
	unsigned u, v, n;

	n = (nrounds + 1) << 1;
	for (u = 0, v = 0; u < n; u++, v += 4) {
		uint64_t x0, x1, x2, x3;

		x0 = x1 = x2 = x3 = comp_skey[u];
		x0 &= (uint64_t)0x1111111111111111;
		x1 &= (uint64_t)0x2222222222222222;
		x2 &= (uint64_t)0x4444444444444444;
		x3 &= (uint64_t)0x8888888888888888;
		x1 >>= 1;
		x2 >>= 2;
		x3 >>= 3;
		skey[v + 0] = (x0 << 4) - x0;
		skey[v + 1] = (x1 << 4) - x1;
		skey[v + 2] = (x2 << 4) - x2;
		skey[v + 3] = (x3 << 4) - x3;
	}
}

static inline void add_round_key(uint64_t *q, const uint64_t *sk) {
	q[0] ^= sk[0];
	q[1] ^= sk[1];
	q[2] ^= sk[2];
	q[3] ^= sk[3];
	q[4] ^= sk[4];
	q[5] ^= sk[5];
	q[6] ^= sk[6];
	q[7] ^= sk[7];
}

static inline void shift_rows(uint64_t *q) {
	int i;

	for (i = 0; i < 8; i++) {
		uint64_t x;

		x = q[i];
		q[i] = (x & (uint64_t)0x000000000000FFFF) | ((x & (uint64_t)0x00000000FFF00000) >> 4) |
		       ((x & (uint64_t)0x00000000000F0000) << 12) | ((x & (uint64_t)0x0000FF0000000000) >> 8) |
		       ((x & (uint64_t)0x000000FF00000000) << 8) | ((x & (uint64_t)0xF000000000000000) >> 12) |
		       ((x & (uint64_t)0x0FFF000000000000) << 4);
	}
}

static inline uint64_t rotr32(uint64_t x) {
	return (x << 32) | (x >> 32);
}

static inline void mix_columns(uint64_t *q) {
	uint64_t q0, q1, q2, q3, q4, q5, q6, q7;
	uint64_t r0, r1, r2, r3, r4, r5, r6, r7;

	q0 = q[0];
	q1 = q[1];
	q2 = q[2];
	q3 = q[3];
	q4 = q[4];
	q5 = q[5];
	q6 = q[6];
	q7 = q[7];
	r0 = (q0 >> 16) | (q0 << 48);
	r1 = (q1 >> 16) | (q1 << 48);
	r2 = (q2 >> 16) | (q2 << 48);
	r3 = (q3 >> 16) | (q3 << 48);
	r4 = (q4 >> 16) | (q4 << 48);
	r5 = (q5 >> 16) | (q5 << 48);
	r6 = (q6 >> 16) | (q6 << 48);
	r7 = (q7 >> 16) | (q7 << 48);

	q[0] = q7 ^ r7 ^ r0 ^ rotr32(q0 ^ r0);
	q[1] = q0 ^ r0 ^ q7 ^ r7 ^ r1 ^ rotr32(q1 ^ r1);
	q[2] = q1 ^ r1 ^ r2 ^ rotr32(q2 ^ r2);
	q[3] = q2 ^ r2 ^ q7 ^ r7 ^ r3 ^ rotr32(q3 ^ r3);
	q[4] = q3 ^ r3 ^ q7 ^ r7 ^ r4 ^ rotr32(q4 ^ r4);
	q[5] = q4 ^ r4 ^ r5 ^ rotr32(q5 ^ r5);
	q[6] = q5 ^ r5 ^ r6 ^ rotr32(q6 ^ r6);
	q[7] = q6 ^ r6 ^ r7 ^ rotr32(q7 ^ r7);
}

static void inc4_be(uint32_t *x) {
	uint32_t t = br_swap32(*x) + 4;
	*x = br_swap32(t);
}

static void aes_ecb4x(unsigned char out[64], const uint32_t ivw[16], const uint64_t *sk_exp, unsigned int nrounds) {
	uint32_t w[16];
	uint64_t q[8];
	unsigned int i;

	memcpy(w, ivw, sizeof(w));
	for (i = 0; i < 4; i++) {
		br_aes_ct64_interleave_in(&q[i], &q[i + 4], w + (i << 2));
	}
	br_aes_ct64_ortho(q);

	add_round_key(q, sk_exp);
	for (i = 1; i < nrounds; i++) {
		br_aes_ct64_bitslice_Sbox(q);
		shift_rows(q);
		mix_columns(q);
		add_round_key(q, sk_exp + (i << 3));
	}
	br_aes_ct64_bitslice_Sbox(q);
	shift_rows(q);
	add_round_key(q, sk_exp + 8 * nrounds);

	br_aes_ct64_ortho(q);
	for (i = 0; i < 4; i++) {
		br_aes_ct64_interleave_out(w + (i << 2), q[i], q[i + 4]);
	}
	br_range_enc32le(out, w, 16);
}

static void aes_ctr4x(unsigned char out[64], uint32_t ivw[16], const uint64_t *sk_exp, unsigned int nrounds) {
	aes_ecb4x(out, ivw, sk_exp, nrounds);

	/* Increase counter for next 4 blocks */
	inc4_be(ivw + 3);
	inc4_be(ivw + 7);
	inc4_be(ivw + 11);
	inc4_be(ivw + 15);
}

static void aes_ecb(unsigned char *out, const unsigned char *in, size_t nblocks, const uint64_t *rkeys, unsigned int nrounds) {
	uint32_t blocks[16];
	unsigned char t[64];

	while (nblocks >= 4) {
		br_range_dec32le(blocks, 16, in);
		aes_ecb4x(out, blocks, rkeys, nrounds);
		nblocks -= 4;
		in += 64;
		out += 64;
	}

	if (nblocks) {
		br_range_dec32le(blocks, nblocks * 4, in);
		aes_ecb4x(t, blocks, rkeys, nrounds);
		memcpy(out, t, nblocks * 16);
	}
}

static inline void aes128_ctr_upd_blks(unsigned char *out, size_t outblks, aes128ctx *ctx) {
	uint32_t ivw[16];
	size_t i;
	uint32_t cc;
	uint8_t *iv = ctx->iv;
	uint32_t blocks = (uint32_t)outblks;
	unsigned int nrounds = 10;

	br_range_dec32le(ivw, 4, iv);

	memcpy(ivw + 4, ivw, 3 * sizeof(uint32_t));
	memcpy(ivw + 8, ivw, 3 * sizeof(uint32_t));
	memcpy(ivw + 12, ivw, 3 * sizeof(uint32_t));
	cc = br_swap32(ivw[3]);
	ivw[7] = br_swap32(cc + 1);
	ivw[11] = br_swap32(cc + 2);
	ivw[15] = br_swap32(cc + 3);

	while (outblks >= 4) {
		aes_ctr4x(out, ivw, ctx->sk_exp, nrounds);
		out += 64;
		outblks -= 4;
	}
	if (outblks > 0) {
		unsigned char tmp[64];
		aes_ctr4x(tmp, ivw, ctx->sk_exp, nrounds);
		for (i = 0; i < outblks * 16; i++) {
			out[i] = tmp[i];
		}
	}
	br_enc32be(&ctx->iv[12], cc + blocks);
}

static inline void aes256_ctr_upd_blks(unsigned char *out, size_t outblks, aes256ctx *ctx) {
	uint32_t ivw[16];
	size_t i;
	uint32_t cc;
	uint8_t *iv = ctx->iv;
	uint32_t blocks = (uint32_t)outblks;
	unsigned int nrounds = 14;

	br_range_dec32le(ivw, 4, iv);

	memcpy(ivw + 4, ivw, 3 * sizeof(uint32_t));
	memcpy(ivw + 8, ivw, 3 * sizeof(uint32_t));
	memcpy(ivw + 12, ivw, 3 * sizeof(uint32_t));
	cc = br_swap32(ivw[3]);
	ivw[7] = br_swap32(cc + 1);
	ivw[11] = br_swap32(cc + 2);
	ivw[15] = br_swap32(cc + 3);

	while (outblks >= 4) {
		aes_ctr4x(out, ivw, ctx->sk_exp, nrounds);
		out += 64;
		outblks -= 4;
	}
	if (outblks > 0) {
		unsigned char tmp[64];
		aes_ctr4x(tmp, ivw, ctx->sk_exp, nrounds);
		for (i = 0; i < outblks * 16; i++) {
			out[i] = tmp[i];
		}
	}
	br_enc32be(&ctx->iv[12], cc + blocks);
}

static void aes_ctr(unsigned char *out, size_t outlen, const unsigned char *iv, const size_t iv_len, const uint64_t *rkeys,
                    unsigned int nrounds) {
	uint32_t ivw[16];
	size_t i;
	uint32_t cc;

	if (iv_len == 12) {
		br_range_dec32le(ivw, 3, iv);
		ivw[3] = 0;
	} else if (iv_len == 16) {
		br_range_dec32le(ivw, 4, iv);
	} else {
		exit(EXIT_FAILURE);
	}
	memcpy(ivw + 4, ivw, 3 * sizeof(uint32_t));
	memcpy(ivw + 8, ivw, 3 * sizeof(uint32_t));
	memcpy(ivw + 12, ivw, 3 * sizeof(uint32_t));
	cc = br_swap32(ivw[3]);
	ivw[7] = br_swap32(cc + 1);
	ivw[11] = br_swap32(cc + 2);
	ivw[15] = br_swap32(cc + 3);

	while (outlen >= 64) {
		aes_ctr4x(out, ivw, rkeys, nrounds);
		out += 64;
		outlen -= 64;
	}
	if (outlen > 0) {
		unsigned char tmp[64];
		aes_ctr4x(tmp, ivw, rkeys, nrounds);
		for (i = 0; i < outlen; i++) {
			out[i] = tmp[i];
		}
	}
}

static void oqs_aes128_load_schedule_c(const uint8_t *key, void **_schedule) {
	*_schedule = malloc(sizeof(aes128ctx));
	OQS_EXIT_IF_NULLPTR(*_schedule, "AES");
	aes128ctx *ctx = (aes128ctx *)*_schedule;
	uint64_t skey[22];
	br_aes_ct64_keysched(skey, key, 16);
	br_aes_ct64_skey_expand(ctx->sk_exp, skey, 10);
}

static void oqs_aes256_load_schedule_c(const uint8_t *key, void **_schedule) {
	*_schedule = malloc(sizeof(aes256ctx));
	OQS_EXIT_IF_NULLPTR(*_schedule, "AES");
	aes256ctx *ctx = (aes256ctx *)*_schedule;
	uint64_t skey[30];
	br_aes_ct64_keysched(skey, key, 32);
	br_aes_ct64_skey_expand(ctx->sk_exp, skey, 14);
}

static void aes_keysched_no_bitslice(uint32_t *skey, const unsigned char *key, unsigned int key_len) {
	unsigned int i, j, k, nk, nkf;
	uint32_t tmp;
	unsigned nrounds = 10 + ((key_len - 16) >> 2);

	nk = (key_len >> 2);
	nkf = ((nrounds + 1) << 2);
	br_range_dec32le(skey, (key_len >> 2), key);
	tmp = skey[(key_len >> 2) - 1];
	for (i = nk, j = 0, k = 0; i < nkf; i++) {
		if (j == 0) {
			tmp = (tmp << 24) | (tmp >> 8);
			tmp = sub_word(tmp) ^ Rcon[k];
		} else if (nk > 6 && j == 4) {
			tmp = sub_word(tmp);
		}
		tmp ^= skey[i - nk];
		skey[i] = tmp;
		if (++j == nk) {
			j = 0;
			k++;
		}
	}
}

static void oqs_aes256_load_schedule_no_bitslice(const uint8_t *key, void **_schedule) {
	*_schedule = malloc(sizeof(aes256ctx_nobitslice));
	assert(*_schedule != NULL);
	uint32_t *schedule = ((aes256ctx_nobitslice *)*_schedule)->sk_exp;
	aes_keysched_no_bitslice(schedule, (const unsigned char *)key, 32);
}

static void oqs_aes256_load_iv_c(const uint8_t *iv, size_t iv_len, void *_schedule) {
	aes256ctx *ctx = _schedule;
	if (iv_len == 12) {
		memcpy(ctx->iv, iv, 12);
		memset(&ctx->iv[12], 0, 4);
	} else if (iv_len == 16) {
		memcpy(ctx->iv, iv, 16);
	} else {
		exit(EXIT_FAILURE);
	}
}

static void oqs_aes256_load_iv_u64_c(uint64_t iv, void *schedule) {
	OQS_EXIT_IF_NULLPTR(schedule, "AES");
	aes256ctx *ctx = (aes256ctx *)schedule;
	ctx->iv[7] = (unsigned char)(iv >> 56);
	ctx->iv[6] = (unsigned char)(iv >> 48);
	ctx->iv[5] = (unsigned char)(iv >> 40);
	ctx->iv[4] = (unsigned char)(iv >> 32);
	ctx->iv[3] = (unsigned char)(iv >> 24);
	ctx->iv[2] = (unsigned char)(iv >> 16);
	ctx->iv[1] = (unsigned char)(iv >> 8);
	ctx->iv[0] = (unsigned char)iv;
	memset(&ctx->iv[8], 0, 8);
}

static void oqs_aes128_load_schedule_no_bitslice(const uint8_t *key, void **_schedule) {
	*_schedule = malloc(44 * sizeof(int));
	assert(*_schedule != NULL);
	uint32_t *schedule = (uint32_t *)*_schedule;
	aes_keysched_no_bitslice(schedule, (const unsigned char *)key, 16);
}

static void oqs_aes128_load_iv_c(const uint8_t *iv, size_t iv_len, void *_schedule) {
	aes128ctx *ctx = _schedule;
	if (iv_len == 12) {
		memcpy(ctx->iv, iv, 12);
		memset(&ctx->iv[12], 0, 4);
	} else if (iv_len == 16) {
		memcpy(ctx->iv, iv, 16);
	} else {
		exit(EXIT_FAILURE);
	}
}

static void oqs_aes128_load_iv_u64_c(uint64_t iv, void *schedule) {
	OQS_EXIT_IF_NULLPTR(schedule, "AES");
	aes128ctx *ctx = (aes128ctx *)schedule;
	ctx->iv[7] = (unsigned char)(iv >> 56);
	ctx->iv[6] = (unsigned char)(iv >> 48);
	ctx->iv[5] = (unsigned char)(iv >> 40);
	ctx->iv[4] = (unsigned char)(iv >> 32);
	ctx->iv[3] = (unsigned char)(iv >> 24);
	ctx->iv[2] = (unsigned char)(iv >> 16);
	ctx->iv[1] = (unsigned char)(iv >> 8);
	ctx->iv[0] = (unsigned char)iv;
	memset(&ctx->iv[8], 0, 8);
}

static void oqs_aes128_ecb_enc_sch_c(const uint8_t *plaintext, const size_t plaintext_len, const void *schedule,
                                     uint8_t *ciphertext) {
	assert(plaintext_len % 16 == 0);
	const aes128ctx *ctx = (const aes128ctx *)schedule;
	aes_ecb(ciphertext, plaintext, plaintext_len / 16, ctx->sk_exp, 10);
}

static void oqs_aes128_ctr_enc_sch_c(const uint8_t *iv, const size_t iv_len, const void *schedule, uint8_t *out,
                                     size_t out_len) {
	const aes128ctx *ctx = (const aes128ctx *)schedule;
	aes_ctr(out, out_len, iv, iv_len, ctx->sk_exp, 10);
}

static void oqs_aes128_ctr_enc_sch_upd_blks_c(void *schedule, uint8_t *out, size_t out_blks) {
	aes128ctx *ctx = (aes128ctx *)schedule;
	aes128_ctr_upd_blks(out, out_blks, ctx);
}

static void oqs_aes256_ecb_enc_sch_c(const uint8_t *plaintext, const size_t plaintext_len, const void *schedule,
                                     uint8_t *ciphertext) {
	assert(plaintext_len % 16 == 0);
	const aes256ctx *ctx = (const aes256ctx *)schedule;
	aes_ecb(ciphertext, plaintext, plaintext_len / 16, ctx->sk_exp, 14);
}

static void oqs_aes256_ctr_enc_sch_c(const uint8_t *iv, const size_t iv_len, const void *schedule, uint8_t *out,
                                     size_t out_len) {
	const aes256ctx *ctx = (const aes256ctx *)schedule;
	aes_ctr(out, out_len, iv, iv_len, ctx->sk_exp, 14);
}

static void oqs_aes256_ctr_enc_sch_upd_blks_c(void *schedule, uint8_t *out, size_t out_blks) {
	aes256ctx *ctx = (aes256ctx *)schedule;
	aes256_ctr_upd_blks(out, out_blks, ctx);
}

static void oqs_aes128_free_schedule_c(void *schedule) {
	if (schedule != NULL) {
		aes128ctx *ctx = (aes128ctx *)schedule;
		OQS_MEM_secure_free(ctx, sizeof(aes128ctx));
	}
}

static void oqs_aes256_free_schedule_c(void *schedule) {
	if (schedule != NULL) {
		aes256ctx *ctx = (aes256ctx *)schedule;
		OQS_MEM_secure_free(ctx, sizeof(aes256ctx));
	}
}

static void oqs_aes256_free_schedule_no_bitslice(void *schedule) {
	if (schedule != NULL) {
		OQS_MEM_secure_free(schedule, sizeof(aes256ctx_nobitslice));
	}
}

static void oqs_aes128_free_schedule_no_bitslice(void *schedule) {
	if (schedule != NULL) {
		OQS_MEM_secure_free(schedule, 44 * sizeof(int));
	}
}

#if __AVX2__

// SPDX-License-Identifier: Public domain
// Based on public domain code by Romain Dolbeau
// http://dolbeau.name/dolbeau/crypto/crypto.html

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
// #include "oqs_common.h"

#include <tmmintrin.h>
#include <wmmintrin.h>

typedef struct {
	__m128i sk_exp[11];
	__m128i iv;
} aes128ctx_;

// From crypto_core/aes128ncrypt/dolbeau/aesenc-int
static inline void aes128ni_setkey_encrypt(const unsigned char *key, __m128i rkeys[11]) {
	__m128i key0 = _mm_loadu_si128((const __m128i *)(key + 0));
	__m128i temp0, temp1, temp4;
	int idx = 0;

	temp0 = key0;

	/* blockshift-based block by Cedric Bourrasset */
#define BLOCK1(IMM)                                \
    temp1 = _mm_aeskeygenassist_si128(temp0, IMM); \
    rkeys[idx++] = temp0;                          \
    temp4 = _mm_slli_si128(temp0, 4);              \
    temp0 = _mm_xor_si128(temp0, temp4);           \
    temp4 = _mm_slli_si128(temp0, 8);              \
    temp0 = _mm_xor_si128(temp0, temp4);           \
    temp1 = _mm_shuffle_epi32(temp1, 0xff);        \
    temp0 = _mm_xor_si128(temp0, temp1)

	BLOCK1(0x01);
	BLOCK1(0x02);
	BLOCK1(0x04);
	BLOCK1(0x08);
	BLOCK1(0x10);
	BLOCK1(0x20);
	BLOCK1(0x40);
	BLOCK1(0x80);
	BLOCK1(0x1b);
	BLOCK1(0x36);
	rkeys[idx++] = temp0;
}

static void oqs_aes128_load_schedule_ni(const uint8_t *key, void **_schedule) {
	*_schedule = malloc(sizeof(aes128ctx_));
	OQS_EXIT_IF_NULLPTR(*_schedule, "AES");
	assert(*_schedule != NULL);
	__m128i *schedule = ((aes128ctx_ *)*_schedule)->sk_exp;
	aes128ni_setkey_encrypt(key, schedule);
}

static void oqs_aes128_load_iv_ni(const uint8_t *iv, size_t iv_len, void *_schedule) {
	aes128ctx_ *ctx = _schedule;
	__m128i idx = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
	if (iv_len == 12) {
		const int32_t *ivi = (const int32_t *)iv;
		ctx->iv = _mm_shuffle_epi8(_mm_set_epi32(0, ivi[2], ivi[1], ivi[0]), idx);
	} else if (iv_len == 16) {
		ctx->iv = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)iv), idx);
	} else {
		exit(EXIT_FAILURE);
	}
}

static void oqs_aes128_load_iv_u64_ni(uint64_t iv, void *_schedule) {
	aes128ctx_ *ctx = _schedule;
	ctx->iv = _mm_loadl_epi64((__m128i *)&iv);
}

static void oqs_aes128_free_schedule_ni(void *schedule) {
	if (schedule != NULL) {
		OQS_MEM_secure_free(schedule, sizeof(aes128ctx_));
	}
}

// From crypto_core/aes128encrypt/dolbeau/aesenc-int
static inline void aes128ni_encrypt(const __m128i rkeys[11], __m128i nv, unsigned char *out) {
	__m128i temp = _mm_xor_si128(nv, rkeys[0]);
	temp = _mm_aesenc_si128(temp, rkeys[1]);
	temp = _mm_aesenc_si128(temp, rkeys[2]);
	temp = _mm_aesenc_si128(temp, rkeys[3]);
	temp = _mm_aesenc_si128(temp, rkeys[4]);
	temp = _mm_aesenc_si128(temp, rkeys[5]);
	temp = _mm_aesenc_si128(temp, rkeys[6]);
	temp = _mm_aesenc_si128(temp, rkeys[7]);
	temp = _mm_aesenc_si128(temp, rkeys[8]);
	temp = _mm_aesenc_si128(temp, rkeys[9]);
	temp = _mm_aesenclast_si128(temp, rkeys[10]);
	_mm_storeu_si128((__m128i *)(out), temp);
}

// 4x interleaved encryption
static inline void aes128ni_encrypt_x4(const __m128i rkeys[11], __m128i n0, __m128i n1, __m128i n2, __m128i n3,
                                       unsigned char *out) {
	__m128i temp0 = _mm_xor_si128(n0, rkeys[0]);
	__m128i temp1 = _mm_xor_si128(n1, rkeys[0]);
	__m128i temp2 = _mm_xor_si128(n2, rkeys[0]);
	__m128i temp3 = _mm_xor_si128(n3, rkeys[0]);

#define AESNENCX4(IDX)                           \
    temp0 = _mm_aesenc_si128(temp0, rkeys[IDX]); \
    temp1 = _mm_aesenc_si128(temp1, rkeys[IDX]); \
    temp2 = _mm_aesenc_si128(temp2, rkeys[IDX]); \
    temp3 = _mm_aesenc_si128(temp3, rkeys[IDX])

	AESNENCX4(1);
	AESNENCX4(2);
	AESNENCX4(3);
	AESNENCX4(4);
	AESNENCX4(5);
	AESNENCX4(6);
	AESNENCX4(7);
	AESNENCX4(8);
	AESNENCX4(9);

	temp0 = _mm_aesenclast_si128(temp0, rkeys[10]);
	temp1 = _mm_aesenclast_si128(temp1, rkeys[10]);
	temp2 = _mm_aesenclast_si128(temp2, rkeys[10]);
	temp3 = _mm_aesenclast_si128(temp3, rkeys[10]);

	_mm_storeu_si128((__m128i *)(out + 0), temp0);
	_mm_storeu_si128((__m128i *)(out + 16), temp1);
	_mm_storeu_si128((__m128i *)(out + 32), temp2);
	_mm_storeu_si128((__m128i *)(out + 48), temp3);
}

static void oqs_aes128_enc_sch_block_ni(const uint8_t *plaintext, const void *_schedule, uint8_t *ciphertext) {
	const __m128i *schedule = ((const aes128ctx_ *)_schedule)->sk_exp;
	aes128ni_encrypt(schedule, _mm_loadu_si128((const __m128i *)plaintext), ciphertext);
}

static void oqs_aes128_ecb_enc_sch_ni(const uint8_t *plaintext, const size_t plaintext_len, const void *schedule,
                                      uint8_t *ciphertext) {
	assert(plaintext_len % 16 == 0);
	for (size_t block = 0; block < plaintext_len / 16; block++) {
		oqs_aes128_enc_sch_block_ni(plaintext + (16 * block), schedule, ciphertext + (16 * block));
	}
}

static void oqs_aes128_ctr_enc_sch_upd_blks_ni(void *schedule, uint8_t *out, size_t out_blks) {
	aes128ctx_ *ctx = (aes128ctx_ *)schedule;
	const __m128i mask = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);

	while (out_blks >= 4) {
		__m128i nv0 = _mm_shuffle_epi8(ctx->iv, mask);
		__m128i nv1 = _mm_shuffle_epi8(_mm_add_epi64(ctx->iv, _mm_set_epi64x(1, 0)), mask);
		__m128i nv2 = _mm_shuffle_epi8(_mm_add_epi64(ctx->iv, _mm_set_epi64x(2, 0)), mask);
		__m128i nv3 = _mm_shuffle_epi8(_mm_add_epi64(ctx->iv, _mm_set_epi64x(3, 0)), mask);
		aes128ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
		ctx->iv = _mm_add_epi64(ctx->iv, _mm_set_epi64x(4, 0));
		out += 64;
		out_blks -= 4;
	}
	while (out_blks >= 1) {
		__m128i nv0 = _mm_shuffle_epi8(ctx->iv, mask);
		aes128ni_encrypt(schedule, nv0, out);
		ctx->iv = _mm_add_epi64(ctx->iv, _mm_set_epi64x(1, 0));
		out += 16;
		out_blks--;
	}
}

static void oqs_aes128_ctr_enc_sch_ni(const uint8_t *iv, const size_t iv_len, const void *schedule, uint8_t *out,
                                      size_t out_len) {
	__m128i block;
	__m128i mask = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
	if (iv_len == 12) {
		const int32_t *ivi = (const int32_t *)iv;
		block = _mm_set_epi32(0, ivi[2], ivi[1], ivi[0]);
	} else if (iv_len == 16) {
		block = _mm_loadu_si128((const __m128i *)iv);
	} else {
		exit(EXIT_FAILURE);
	}

	while (out_len >= 64) {
		__m128i nv0 = block;
		__m128i nv1 = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)), mask);
		__m128i nv2 = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(2, 0)), mask);
		__m128i nv3 = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(3, 0)), mask);
		aes128ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
		block = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(4, 0)), mask);
		out += 64;
		out_len -= 64;
	}
	while (out_len >= 16) {
		aes128ni_encrypt(schedule, block, out);
		out += 16;
		out_len -= 16;
		block = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)), mask);
	}
	if (out_len > 0) {
		uint8_t tmp[16];
		aes128ni_encrypt(schedule, block, tmp);
		memcpy(out, tmp, out_len);
	}
}

// SPDX-License-Identifier: Public domain
// Based on public domain code by Romain Dolbeau
// http://dolbeau.name/dolbeau/crypto/crypto.html

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
// #include "oqs_common.h"

#include <tmmintrin.h>
#include <wmmintrin.h>

#define AES_BLOCKBYTES 16

typedef struct {
	__m128i sk_exp[15];
	__m128i iv;
} aes256ctx_;

#define BE_TO_UINT32(n)                                                                                    \
    (uint32_t)((((uint8_t *)&(n))[0] << 24) | (((uint8_t *)&(n))[1] << 16) | (((uint8_t *)&(n))[2] << 8) | \
               (((uint8_t *)&(n))[3] << 0))

// From crypto_core/aes256encrypt/dolbeau/aesenc-int
static inline void aes256ni_setkey_encrypt(const unsigned char *key, __m128i rkeys[15]) {
	__m128i key0 = _mm_loadu_si128((const __m128i *)(key + 0));
	__m128i key1 = _mm_loadu_si128((const __m128i *)(key + 16));
	__m128i temp0, temp1, temp2, temp4;
	int idx = 0;

	rkeys[idx++] = key0;
	temp0 = key0;
	temp2 = key1;

	/* blockshift-based block by Cedric Bourrasset & Romain Dolbeau */
#undef BLOCK1
#define BLOCK1(IMM)                                \
    temp1 = _mm_aeskeygenassist_si128(temp2, IMM); \
    rkeys[idx++] = temp2;                          \
    temp4 = _mm_slli_si128(temp0, 4);              \
    temp0 = _mm_xor_si128(temp0, temp4);           \
    temp4 = _mm_slli_si128(temp0, 8);              \
    temp0 = _mm_xor_si128(temp0, temp4);           \
    temp1 = _mm_shuffle_epi32(temp1, 0xff);        \
    temp0 = _mm_xor_si128(temp0, temp1)

#define BLOCK2(IMM)                                \
    temp1 = _mm_aeskeygenassist_si128(temp0, IMM); \
    rkeys[idx++] = temp0;                          \
    temp4 = _mm_slli_si128(temp2, 4);              \
    temp2 = _mm_xor_si128(temp2, temp4);           \
    temp4 = _mm_slli_si128(temp2, 8);              \
    temp2 = _mm_xor_si128(temp2, temp4);           \
    temp1 = _mm_shuffle_epi32(temp1, 0xaa);        \
    temp2 = _mm_xor_si128(temp2, temp1)

	BLOCK1(0x01);
	BLOCK2(0x01);

	BLOCK1(0x02);
	BLOCK2(0x02);

	BLOCK1(0x04);
	BLOCK2(0x04);

	BLOCK1(0x08);
	BLOCK2(0x08);

	BLOCK1(0x10);
	BLOCK2(0x10);

	BLOCK1(0x20);
	BLOCK2(0x20);

	BLOCK1(0x40);
	rkeys[idx++] = temp0;
}

static void oqs_aes256_load_schedule_ni(const uint8_t *key, void **_schedule) {
	*_schedule = malloc(sizeof(aes256ctx_));
	OQS_EXIT_IF_NULLPTR(*_schedule, "AES");
	assert(*_schedule != NULL);
	__m128i *schedule = ((aes256ctx_ *)*_schedule)->sk_exp;
	aes256ni_setkey_encrypt(key, schedule);
}

static void oqs_aes256_load_iv_ni(const uint8_t *iv, size_t iv_len, void *_schedule) {
	aes256ctx_ *ctx = _schedule;
	__m128i idx = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
	if (iv_len == 12) {
		const int32_t *ivi = (const int32_t *)iv;
		ctx->iv = _mm_shuffle_epi8(_mm_set_epi32(0, ivi[2], ivi[1], ivi[0]), idx);
	} else if (iv_len == 16) {
		ctx->iv = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i *)iv), idx);
	} else {
		exit(EXIT_FAILURE);
	}
}

static void oqs_aes256_load_iv_u64_ni(uint64_t iv, void *_schedule) {
	aes256ctx_ *ctx = _schedule;
	ctx->iv = _mm_loadl_epi64((__m128i *)&iv);
}

static void oqs_aes256_free_schedule_ni(void *schedule) {
	if (schedule != NULL) {
		OQS_MEM_secure_free(schedule, sizeof(aes256ctx_));
	}
}

// Modified from crypto_core/aes256encrypt/dolbeau/aesenc-int
static inline void aes256ni_encrypt(const __m128i rkeys[15], __m128i nv, unsigned char *out) {
	__m128i temp = _mm_xor_si128(nv, rkeys[0]);
	temp = _mm_aesenc_si128(temp, rkeys[1]);
	temp = _mm_aesenc_si128(temp, rkeys[2]);
	temp = _mm_aesenc_si128(temp, rkeys[3]);
	temp = _mm_aesenc_si128(temp, rkeys[4]);
	temp = _mm_aesenc_si128(temp, rkeys[5]);
	temp = _mm_aesenc_si128(temp, rkeys[6]);
	temp = _mm_aesenc_si128(temp, rkeys[7]);
	temp = _mm_aesenc_si128(temp, rkeys[8]);
	temp = _mm_aesenc_si128(temp, rkeys[9]);
	temp = _mm_aesenc_si128(temp, rkeys[10]);
	temp = _mm_aesenc_si128(temp, rkeys[11]);
	temp = _mm_aesenc_si128(temp, rkeys[12]);
	temp = _mm_aesenc_si128(temp, rkeys[13]);
	temp = _mm_aesenclast_si128(temp, rkeys[14]);
	_mm_storeu_si128((__m128i *)(out), temp);
}

// 4x interleaved encryption
static inline void aes256ni_encrypt_x4(const __m128i rkeys[15], __m128i n0, __m128i n1, __m128i n2, __m128i n3,
                                       unsigned char *out) {
	__m128i temp0 = _mm_xor_si128(n0, rkeys[0]);
	__m128i temp1 = _mm_xor_si128(n1, rkeys[0]);
	__m128i temp2 = _mm_xor_si128(n2, rkeys[0]);
	__m128i temp3 = _mm_xor_si128(n3, rkeys[0]);

#define AESNENCX4(IDX)                           \
    temp0 = _mm_aesenc_si128(temp0, rkeys[IDX]); \
    temp1 = _mm_aesenc_si128(temp1, rkeys[IDX]); \
    temp2 = _mm_aesenc_si128(temp2, rkeys[IDX]); \
    temp3 = _mm_aesenc_si128(temp3, rkeys[IDX])

	AESNENCX4(1);
	AESNENCX4(2);
	AESNENCX4(3);
	AESNENCX4(4);
	AESNENCX4(5);
	AESNENCX4(6);
	AESNENCX4(7);
	AESNENCX4(8);
	AESNENCX4(9);
	AESNENCX4(10);
	AESNENCX4(11);
	AESNENCX4(12);
	AESNENCX4(13);

	temp0 = _mm_aesenclast_si128(temp0, rkeys[14]);
	temp1 = _mm_aesenclast_si128(temp1, rkeys[14]);
	temp2 = _mm_aesenclast_si128(temp2, rkeys[14]);
	temp3 = _mm_aesenclast_si128(temp3, rkeys[14]);

	_mm_storeu_si128((__m128i *)(out + 0), temp0);
	_mm_storeu_si128((__m128i *)(out + 16), temp1);
	_mm_storeu_si128((__m128i *)(out + 32), temp2);
	_mm_storeu_si128((__m128i *)(out + 48), temp3);
}

static void oqs_aes256_enc_sch_block_ni(const uint8_t *plaintext, const void *_schedule, uint8_t *ciphertext) {
	const __m128i *schedule = ((const aes256ctx_ *)_schedule)->sk_exp;
	aes256ni_encrypt(schedule, _mm_loadu_si128((const __m128i *)plaintext), ciphertext);
}

static void oqs_aes256_ecb_enc_sch_ni(const uint8_t *plaintext, const size_t plaintext_len, const void *schedule,
                                      uint8_t *ciphertext) {
	assert(plaintext_len % 16 == 0);
	for (size_t block = 0; block < plaintext_len / 16; block++) {
		oqs_aes256_enc_sch_block_ni(plaintext + (16 * block), schedule, ciphertext + (16 * block));
	}
}

static void oqs_aes256_ctr_enc_sch_upd_blks_ni(void *schedule, uint8_t *out, size_t out_blks) {
	aes256ctx_ *ctx = (aes256ctx_ *)schedule;
	const __m128i mask = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);

	while (out_blks >= 4) {
		__m128i nv0 = _mm_shuffle_epi8(ctx->iv, mask);
		__m128i nv1 = _mm_shuffle_epi8(_mm_add_epi64(ctx->iv, _mm_set_epi64x(1, 0)), mask);
		__m128i nv2 = _mm_shuffle_epi8(_mm_add_epi64(ctx->iv, _mm_set_epi64x(2, 0)), mask);
		__m128i nv3 = _mm_shuffle_epi8(_mm_add_epi64(ctx->iv, _mm_set_epi64x(3, 0)), mask);
		aes256ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
		ctx->iv = _mm_add_epi64(ctx->iv, _mm_set_epi64x(4, 0));
		out += 64;
		out_blks -= 4;
	}
	while (out_blks >= 1) {
		__m128i nv0 = _mm_shuffle_epi8(ctx->iv, mask);
		aes256ni_encrypt(schedule, nv0, out);
		ctx->iv = _mm_add_epi64(ctx->iv, _mm_set_epi64x(1, 0));
		out += 16;
		out_blks--;
	}
}

static void oqs_aes256_ctr_enc_sch_ni(const uint8_t *iv, const size_t iv_len, const void *schedule, uint8_t *out,
                                      size_t out_len) {
	__m128i block;
	__m128i mask = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);
	if (iv_len == 12) {
		const int32_t *ivi = (const int32_t *)iv;
		block = _mm_set_epi32(0, ivi[2], ivi[1], ivi[0]);
	} else if (iv_len == 16) {
		block = _mm_loadu_si128((const __m128i *)iv);
	} else {
		exit(EXIT_FAILURE);
	}

	while (out_len >= 64) {
		__m128i nv0 = block;
		__m128i nv1 = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)), mask);
		__m128i nv2 = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(2, 0)), mask);
		__m128i nv3 = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(3, 0)), mask);
		aes256ni_encrypt_x4(schedule, nv0, nv1, nv2, nv3, out);
		block = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(4, 0)), mask);
		out += 64;
		out_len -= 64;
	}
	while (out_len >= 16) {
		aes256ni_encrypt(schedule, block, out);
		out += 16;
		out_len -= 16;
		block = _mm_shuffle_epi8(_mm_add_epi64(_mm_shuffle_epi8(block, mask), _mm_set_epi64x(1, 0)), mask);
	}
	if (out_len > 0) {
		uint8_t tmp[16];
		aes256ni_encrypt(schedule, block, tmp);
		memcpy(out, tmp, out_len);
	}
}

#endif

/**
 * Add missing function needed by the oqs derived aes
 *
 * Adapted from
 *
 * https://github.com/open-quantum-safe/liboqs
 * commit 3488f0a598c64b730ee2e2a4acb38e1a51797c99
 *
 * MIT license
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#if __AVX2__

int AES128_CTR(unsigned char *output, size_t outputByteLen, const unsigned char *input, size_t inputByteLen) {
	assert(inputByteLen > 0);
	const uint8_t iv[16] = {0};
	void *schedule = NULL;
	oqs_aes128_load_schedule_ni(input, &schedule);
	oqs_aes128_ctr_enc_sch_ni(iv, 16, schedule, output, outputByteLen);
	oqs_aes128_free_schedule_ni(schedule);
	return (int)outputByteLen;
}

void AES256_ECB(const unsigned char *key, const uint8_t *input, unsigned char *output) {
	void *schedule = NULL;
	oqs_aes256_load_schedule_ni(key, &schedule);
	oqs_aes256_ecb_enc_sch_ni(input, 16, schedule, output);
	oqs_aes256_free_schedule_ni(schedule);
}

#else

int AES128_CTR(unsigned char *output, size_t outputByteLen, const unsigned char *input, size_t inputByteLen) {
	assert(inputByteLen > 0);
	const uint8_t iv[16] = {0};
	void *schedule = NULL;
	oqs_aes128_load_schedule_c(input, &schedule);
	oqs_aes128_ctr_enc_sch_c(iv, 16, schedule, output, outputByteLen);
	oqs_aes128_free_schedule_c(schedule);
	return (int)outputByteLen;
}

void AES256_ECB(const unsigned char *key, const uint8_t *input, unsigned char *output) {
	void *schedule = NULL;
	oqs_aes256_load_schedule_c(key, &schedule);
	oqs_aes256_ecb_enc_sch_c(input, 16, schedule, output);
	oqs_aes256_free_schedule_c(schedule);
}

#endif
#endif
