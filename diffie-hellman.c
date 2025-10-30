#include <stdio.h>
#include <time.h>
#include <gmp.h>

// trial division to factor n into distinct prime factors (sufficient for primitive-root test)
void factor_distinct(mpz_t n, mpz_t factors[], size_t *k) {
    mpz_t d, tmp;
    mpz_inits(d, tmp, NULL);
    *k = 0;

    // factor out 2
    if (mpz_divisible_ui_p(n, 2)) {
        mpz_init_set_ui(factors[(*k)++], 2);
        while (mpz_divisible_ui_p(n, 2)) mpz_divexact_ui(n, n, 2);
    }
    // odd trial division
    for (unsigned long i = 3; mpz_cmp_ui(n, 1) > 0; i += 2) {
        if (mpz_divisible_ui_p(n, i)) {
            mpz_init_set_ui(factors[(*k)++], i);
            while (mpz_divisible_ui_p(n, i)) mpz_divexact_ui(n, n, i);
        }
        // stop when i*i > n
        mpz_set_ui(tmp, i);
        mpz_mul(tmp, tmp, tmp);
        if (mpz_cmp(tmp, n) > 0) break;
    }
    if (mpz_cmp_ui(n, 1) > 0) { // leftover prime
        mpz_init_set(factors[(*k)++], n);
    }
    mpz_clears(d, tmp, NULL);
}

// return 1 if g is a primitive root mod p (p prime) given factors of p-1
int is_generator(const mpz_t g, const mpz_t p, mpz_t factors[], size_t k) {
    mpz_t p_minus_1, exp, t;
    mpz_inits(p_minus_1, exp, t, NULL);
    mpz_sub_ui(p_minus_1, p, 1);

    for (size_t i = 0; i < k; ++i) {
        mpz_divexact(exp, p_minus_1, factors[i]);
        mpz_powm(t, g, exp, p);
        if (mpz_cmp_ui(t, 1) == 0) { // not a generator
            mpz_clears(p_minus_1, exp, t, NULL);
            return 0;
        }
    }
    mpz_clears(p_minus_1, exp, t, NULL);
    return 1;
}

int main(void) {
    // i) Choose a ≥30-digit prime q = P
    mpz_t P; mpz_init(P);
    // Example 30-digit prime (replace with your chosen prime if needed)
    mpz_set_str(P, "982451653173961852241334935997", 10);

    // ii) Find primitive root α immediately greater than last 2 digits of your UMBC ID.
    // Put your two-digit threshold here:
    unsigned long threshold = 15; // e.g., if your last two digits are 15
    mpz_t Pm1; mpz_init(Pm1); mpz_sub_ui(Pm1, P, 1);

    // factor P-1
    mpz_t tmp; mpz_init_set(tmp, Pm1);
    mpz_t factors[256]; size_t k = 0;
    factor_distinct(tmp, factors, &k);

    // search for generator α >= threshold+1
    mpz_t alpha; mpz_init(alpha);
    clock_t t0 = clock();
    for (unsigned long cand = threshold + 1;; ++cand) {
        mpz_set_ui(alpha, cand);
        if (is_generator(alpha, P, factors, k)) break;
    }
    clock_t t1 = clock();
    double seconds = (double)(t1 - t0) / CLOCKS_PER_SEC;

    // iii) Choose private keys XA, XB (> last 5 digits of your UMBC ID)
    mpz_t XA, XB; mpz_inits(XA, XB, NULL);
    mpz_set_ui(XA, 51015); // replace with your values
    mpz_set_ui(XB, 51016);

    // iv) Compute YA = α^XA mod P, YB = α^XB mod P
    mpz_t YA, YB; mpz_inits(YA, YB, NULL);
    mpz_powm(YA, alpha, XA, P);
    mpz_powm(YB, alpha, XB, P);

    // v) Shared key SAB
    mpz_t SA, SB; mpz_inits(SA, SB, NULL);
    mpz_powm(SA, YB, XA, P);
    mpz_powm(SB, YA, XB, P);

    gmp_printf("P (prime)  = %Zd\n", P);
    gmp_printf("alpha (g)  = %Zd\n", alpha);
    printf("Primitive root search time: %.6f s\n", seconds);
    gmp_printf("XA         = %Zd\n", XA);
    gmp_printf("XB         = %Zd\n", XB);
    gmp_printf("YA         = %Zd\n", YA);
    gmp_printf("YB         = %Zd\n", YB);
    gmp_printf("S_A        = %Zd\n", SA);
    gmp_printf("S_B        = %Zd\n", SB);
    printf("Keys match? %s\n", mpz_cmp(SA, SB) == 0 ? "YES" : "NO");

    // cleanup
    for (size_t i = 0; i < k; ++i) mpz_clear(factors[i]);
    mpz_clears(P, Pm1, tmp, alpha, XA, XB, YA, YB, SA, SB, NULL);
    return 0;
}
