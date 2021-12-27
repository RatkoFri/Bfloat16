#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <inttypes.h>


#define FLOAT_MANT_BITS    (23)
#define FLOAT_EXPO_BITS    (8)
#define FLOAT_EXPO_BIAS    (127)
#define FLOAT_MANT_MASK    (~((~0u) << (FLOAT_MANT_BITS+1))) /* incl. integer bit */
#define EXPO_ADJUST        (1)   /* adjustment for performance reasons */
#define MIN_NORM_EXPO      (1)   /* minimum biased exponent of normals */
#define MAX_NORM_EXPO      (254) /* maximum biased exponent of normals */
#define INF_EXPO           (255) /* biased exponent of infinities */
#define EXPO_MASK          (~((~0u) << FLOAT_EXPO_BITS))
#define FLOAT_SIGN_MASK    (0x80000000u)
#define FLOAT_IMPLICIT_BIT (1 << FLOAT_MANT_BITS)
#define RND_BIT_SHIFT      (31)
#define RND_BIT_MASK       (1u << RND_BIT_SHIFT)
#define FLOAT_INFINITY     (0x7f800000)
#define FLOAT_INDEFINITE   (0xffc00000u)
#define MANT_LSB           (0x00000001)
#define FLOAT_QNAN_BIT     (0x00400000)
#define MAX_SHIFT          (FLOAT_MANT_BITS + 2)


void float2bfloat(const float src, float& dst) {
	const uint16_t* p = reinterpret_cast<const uint16_t*>(&src);
	uint16_t* q = reinterpret_cast<uint16_t*>(&dst);
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  q[0] = p[0];
  q[1] = 0;
#else
	q[0] = 0;
	q[1] = p[1];
#endif
}

uint32_t fp32_mul_core (uint32_t a, uint32_t b)
{
    uint64_t prod;
    uint32_t expoa, expob, manta, mantb, shift;
    uint32_t r, signr, expor, mantr_hi, mantr_lo;

    /* split arguments into sign, exponent, significand */
    expoa = ((a >> FLOAT_MANT_BITS) & EXPO_MASK) - EXPO_ADJUST;
    expob = ((b >> FLOAT_MANT_BITS) & EXPO_MASK) - EXPO_ADJUST;
    manta = (a | FLOAT_IMPLICIT_BIT) & FLOAT_MANT_MASK;
    mantb = (b | FLOAT_IMPLICIT_BIT) & FLOAT_MANT_MASK;
    /* result sign bit: XOR sign argument signs */
    signr = (a ^ b) & FLOAT_SIGN_MASK;
    if ((expoa >= (MAX_NORM_EXPO - EXPO_ADJUST)) || /* at least one argument is special */
        (expob >= (MAX_NORM_EXPO - EXPO_ADJUST))) { 
        if ((a & ~FLOAT_SIGN_MASK) > FLOAT_INFINITY) { /* a is NaN */
            /* return quietened NaN */
            return a | FLOAT_QNAN_BIT;
        }
        if ((b & ~FLOAT_SIGN_MASK) > FLOAT_INFINITY) { /* b is NaN */
            /* return quietened NaN */
            return b | FLOAT_QNAN_BIT;
        }
        if ((a & ~FLOAT_SIGN_MASK) == 0) { /* a is zero */
            /* return NaN if b is infinity, else zero */
            return (expob != (INF_EXPO - EXPO_ADJUST)) ? signr : FLOAT_INDEFINITE;
        }
        if ((b & ~FLOAT_SIGN_MASK) == 0) { /* b is zero */
            /* return NaN if a is infinity, else zero */
            return (expoa != (INF_EXPO - EXPO_ADJUST)) ? signr : FLOAT_INDEFINITE;
        }
        if (((a & ~FLOAT_SIGN_MASK) == FLOAT_INFINITY) || /* a or b infinity */
            ((b & ~FLOAT_SIGN_MASK) == FLOAT_INFINITY)) {
            return signr | FLOAT_INFINITY;
        }
        if ((int32_t)expoa < (MIN_NORM_EXPO - EXPO_ADJUST)) { /* a is subnormal */
            /* normalize significand of a */
            manta = a & FLOAT_MANT_MASK;
            expoa++;
            do {
                manta = 2 * manta;
                expoa--;
            } while (manta < FLOAT_IMPLICIT_BIT);
        } else if ((int32_t)expob < (MIN_NORM_EXPO - EXPO_ADJUST)) { /* b is subnormal */
            /* normalize significand of b */
            mantb = b & FLOAT_MANT_MASK;
            expob++;
            do {
                mantb = 2 * mantb;
                expob--;
            } while (mantb < FLOAT_IMPLICIT_BIT);
        }
    }
    /* result exponent: add argument exponents and adjust for biasing */
    expor = expoa + expob - FLOAT_EXPO_BIAS + 2 * EXPO_ADJUST;
    mantb = mantb ; /* preshift to align result signficand */
    /* result significand: multiply argument signficands */
    prod = (uint64_t)manta * mantb;
    prod = prod << FLOAT_EXPO_BITS;
    mantr_hi = (uint32_t)(prod >> 32);
    mantr_lo = (uint32_t)(prod >>  0);
    /* normalize significand */
    if (mantr_hi < FLOAT_IMPLICIT_BIT) {
        mantr_hi = (mantr_hi << 1) | (mantr_lo >> (32 - 1));
        mantr_lo = (mantr_lo << 1);
        expor--;
    }
    if (expor <= (MAX_NORM_EXPO - EXPO_ADJUST)) { /* normal, may overflow to infinity during rounding */
        /* combine biased exponent, sign and signficand */
        r = (expor << FLOAT_MANT_BITS) + signr + mantr_hi;
        /* round result to nearest or even; overflow to infinity possible */
        r = r + ((mantr_lo == RND_BIT_MASK) ? (mantr_hi & MANT_LSB) : (mantr_lo >> RND_BIT_SHIFT));
    } else if ((int32_t)expor > (MAX_NORM_EXPO - EXPO_ADJUST)) { /* overflow */
        /* return infinity */
        r = signr | FLOAT_INFINITY;
    } else { /* underflow */
        /* return zero, normal, or smallest subnormal */
        shift = 0 - expor;
        if (shift > MAX_SHIFT) shift = MAX_SHIFT;
        /* denormalize significand */
        mantr_lo = mantr_hi << (32 - shift) | (mantr_lo ? 1 : 0);
        mantr_hi = mantr_hi >> shift;
        /* combine sign and signficand; biased exponent known to be zero */
        r = mantr_hi + signr;
        /* round result to nearest or even */
        r = r + ((mantr_lo == RND_BIT_MASK) ? (mantr_hi & MANT_LSB) : (mantr_lo >> RND_BIT_SHIFT));
    }
    return r;
}

uint32_t float_as_uint (float a)
{
    uint32_t r;
    memcpy (&r, &a, sizeof r);
    return r;
}

float uint_as_float (uint32_t a)
{
    float r;
    memcpy (&r, &a, sizeof r);
    return r;
}

float fp32_mul (float a, float b)
{
    return uint_as_float (fp32_mul_core (float_as_uint (a), float_as_uint (b)));
}

/*
Strategija: 

1. Konverzija iz Floating pointa v BFloat16 - f-ja float2bfloat
    Slabosti: Groba konverzija v BFloat16, samo daje stran spodnjih 16 bitov 
2. Mnozenje Bfloat16 stevil s pomocjo single precision mnozilnika. 
    Prednosti: Super koda. Dragac pa pri mnozenju ne smemo pozabiti da so spodnjih 16 bitov vedno enaki 0
    Slabosti: Overkill za CNNove, prevec se ukvarja z nedovoljenim slucaji da je noro. Za simulacijsko verzijo je OK, 
    za CNN pa ne vem. 

*/
#define NUM_RAND 10000
// F-ja za generiranje floating point števil
float RandomFloat(float min, float max){
   return ((max - min) * ((float)rand() / RAND_MAX)) + min;
}

int main (void)
{
   	float a = 2.5;
	float b = 3.1245;
	float tmp_b,tmp_a;

    FILE *log;
	char p;

	float2bfloat(b,tmp_b);
	float2bfloat(a,tmp_a);

	printf("True product %f\n", a*b);
	printf("Reduced precision product %f\n", tmp_a*tmp_b);
	printf("Reduced precision product %f\n", fp32_mul(tmp_a,tmp_b));
    float p_exact, p_approx;

    float RE = 0;
    int i;
    for(i=0;i<RAND_MAX;i++){
        float2bfloat(RandomFloat(-10000,10000),tmp_a);
        float2bfloat(RandomFloat(-10000,10000),tmp_b);
        p_exact = tmp_a * tmp_b;
        p_approx = fp32_mul(tmp_a,tmp_b);
        if(p_exact == 0){
            RE += 0;
        }
        else{
            RE += abs(p_exact - p_approx)/abs(p_exact);
        }
    }

    printf("Relative error %f\n ", RE/RAND_MAX/RAND_MAX);

    log = fopen("log/EX_log.txt", "a");

		if (log != NULL)
		{	
			
			char str0[] = "--------------------------------- \n \n";
			char str1[80];

			fputs(str0, log);
			sprintf(str1, "Exact floating point  \n" );
			fputs(str1, log);

			sprintf(str1, "Average RE %2.10g %% \n", (RE/RAND_MAX/RAND_MAX));
			fputs(str1, log);

			fclose(log);
		}

	return 0; 
}