// diffie_fast.c
// Build: gcc -O2 diffie_fast.c -o diffie_fast -lgmp
// Run  : ./diffie_fast
//
// What it does (fast path only):
// - Generates a safe prime P = 2*r + 1 with > 40 digits (default ~170 bits ≈ 51 digits).
// - Factors of P-1 are just {2, r}; uses that to test primitive roots quickly.
// - Finds a primitive root alpha (≥ 100 to satisfy "≥ 3 digits").
// - Runs DH with sample private exponents (customizable).
// - Times the primitive-root search.
//
// Notes:
// - If you *must* use a specific P and it's *also* safe (i.e., (P-1)/2 is prime),
//   set USE_HARDCODED_P = 1 and fill HARDCODED_P below. Otherwise, keep generation on.
//
// Optional tweaks near the top:
//   DIGITS_MIN   : minimum decimal digits for P (must be ≥ 41)
//   XA_UI / XB_UI: private exponents (must exceed your assignment's thresholds)

#include <stdio.h>
#include <time.h>
#include <gmp.h>

// -------------------- Config --------------------
#define USE_HARDCODED_P 0
// If USE_HARDCODED_P == 1, put a SAFE PRIME here (P = 2r+1, with r prime):
#define HARDCODED_P "0"

// Minimum digits for P (assignment requires > 40). 51 ≈ 170 bits.
static const unsigned DIGITS_MIN = 51;

// Private exponents (edit as needed; > your ID's last 5 digits, per assignment)
static const unsigned long XA_UI = 51015;
static const unsigned long XB_UI = 51016;

// Generator search starts at least from 100 (≥ 3 digits)
static const unsigned long GEN_START_MIN = 100;

// Miller-Rabin reps
static const int PRP_REPS = 30;

// ------------------------------------------------

// Convert approximate decimal digits to bits: digits * log2(10) ≈ digits * 3.32193
static unsigned digits_to_bits(unsigned digits) {
    double bits_d = digits * 3.3219280948873626;
    unsigned bits = (unsigned)(bits_d + 0.5);
    if (bits < 3) bits = 3;
    return bits;
}

// Generate a safe prime P=2r+1 with at least 'digits' decimal digits
static void gen_safe_prime(mpz_t P, mpz_t r, unsigned digits) {
    gmp_randstate_t st;
    gmp_randinit_default(st);
    gmp_randseed_ui(st, (unsigned)time(NULL));

    unsigned bits = digits_to_bits(digits);
    if (bits < 130) bits = 130; // keep it reasonably large

    mpz_t candidate_r;
    mpz_init(candidate_r);

    // Loop until both r and P=2r+1 are probably prime
    for (;;) {
        mpz_urandomb(candidate_r, st, bits - 1);   // r has ~bits-1 bits
        mpz_setbit(candidate_r, bits - 2);         // ensure high bit set for size
        mpz_nextprime(candidate_r, candidate_r);   // r = next prime

        // P = 2r + 1
        mpz_mul_ui(P, candidate_r, 2);
        mpz_add_ui(P, P, 1);

        if (mpz_probab_prime_p(P, PRP_REPS)) {
            // Ensure P has at least 'digits' decimal digits
            if (mpz_sizeinbase(P, 10) >= digits) break;
        }
    }
    mpz_set(r, candidate_r);
    mpz_clear(candidate_r);
    gmp_randclear(st);
}

// Return 1 if g is a primitive root modulo safe prime P (P=2r+1) using factors {2, r}
static int is_generator_safe_prime(const mpz_t g, const mpz_t P, const mpz_t r) {
    // For a primitive root modulo safe prime P, we need:
    // g^((P-1)/2) != 1 mod P   and   g^((P-1)/r) != 1 mod P
    mpz_t exp, t, Pm1;
    mpz_inits(exp, t, Pm1, NULL);
    mpz_sub_ui(Pm1, P, 1);

    // Check factor 2
    mpz_divexact_ui(exp, Pm1, 2);
    mpz_powm(t, g, exp, P);
    if (mpz_cmp_ui(t, 1) == 0) { mpz_clears(exp, t, Pm1, NULL); return 0; }

    // Check factor r
    mpz_divexact(exp, Pm1, r);
    mpz_powm(t, g, exp, P);
    if (mpz_cmp_ui(t, 1) == 0) { mpz_clears(exp, t, Pm1, NULL); return 0; }

    mpz_clears(exp, t, Pm1, NULL);
    return 1;
}

int main(void) {
    mpz_t P, r; mpz_inits(P, r, NULL);

#if USE_HARDCODED_P
    // Use a hardcoded SAFE prime (P=2r+1). Ensure HARDCODED_P is safe!
    if (mpz_set_str(P, HARDCODED_P, 10) != 0) {
        fprintf(stderr, "Invalid HARDCODED_P.\n");
        mpz_clears(P, r, NULL);
        return 1;
    }
    if (!mpz_probab_prime_p(P, PRP_REPS)) {
        fprintf(stderr, "HARDCODED_P is not prime.\n");
        mpz_clears(P, r, NULL);
        return 1;
    }
    // Compute r = (P-1)/2 and check it's prime
    mpz_sub_ui(r, P, 1);
    if (!mpz_divisible_ui_p(r, 2)) {
        fprintf(stderr, "HARDCODED_P is not of the form 2r+1.\n");
        mpz_clears(P, r, NULL);
        return 1;
    }
    mpz_divexact_ui(r, r, 2);
    if (!mpz_probab_prime_p(r, PRP_REPS)) {
        fprintf(stderr, "(P-1)/2 is not prime; P is not safe.\n");
        mpz_clears(P, r, NULL);
        return 1;
    }
#else
    // Generate a safe prime with >= DIGITS_MIN decimal digits
    gen_safe_prime(P, r, DIGITS_MIN);
#endif

    // Primitive root search (start at ≥ GEN_START_MIN)
    mpz_t alpha; mpz_init(alpha);
    clock_t t0 = clock();
    unsigned long start = GEN_START_MIN;
    for (unsigned long cand = start;; ++cand) {
        mpz_set_ui(alpha, cand);
        if (is_generator_safe_prime(alpha, P, r)) break;
    }
    clock_t t1 = clock();
    double seconds = (double)(t1 - t0) / CLOCKS_PER_SEC;

    // Private exponents
    mpz_t XA, XB; mpz_inits(XA, XB, NULL);
    mpz_set_ui(XA, XA_UI);
    mpz_set_ui(XB, XB_UI);

    // Public keys
    mpz_t YA, YB; mpz_inits(YA, YB, NULL);
    mpz_powm(YA, alpha, XA, P);
    mpz_powm(YB, alpha, XB, P);

    // Shared secrets
    mpz_t SA, SB; mpz_inits(SA, SB, NULL);
    mpz_powm(SA, YB, XA, P);
    mpz_powm(SB, YA, XB, P);

    // Output
    gmp_printf("P (prime, %lu digits) = %Zd\n", mpz_sizeinbase(P, 10), P);
    gmp_printf("r ( (P-1)/2, prime ) = %Zd\n", r);
    gmp_printf("alpha (generator)     = %Zd\n", alpha);
    printf("Primitive root search time: %.6f s\n", seconds);
    gmp_printf("XA = %Zd\n", XA);
    gmp_printf("XB = %Zd\n", XB);
    gmp_printf("YA = %Zd\n", YA);
    gmp_printf("YB = %Zd\n", YB);
    gmp_printf("S_A = %Zd\n", SA);
    gmp_printf("S_B = %Zd\n", SB);
    printf("Keys match? %s\n", (mpz_cmp(SA, SB) == 0) ? "YES" : "NO");

    // Cleanup
    mpz_clears(P, r, alpha, XA, XB, YA, YB, SA, SB, NULL);
    return 0;
}
