// https://www.geeksforgeeks.org/computer-networks/rsa-algorithm-cryptography/
// https://www.geeksforgeeks.org/dsa/euclidean-algorithms-basic-and-extended/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

// Function for extended Euclidean Algorithm 
int gcdExtended(int a, int b, int *x, int *y) {

    // Base Case 
    if (a == 0) { 
        *x = 0; 
        *y = 1; 
        return b; 
    } 

    int x1, y1; 
    int gcd = gcdExtended(b % a, a, &x1, &y1); 

    // Update x and y using results of 
    // recursive call 
    *x = y1 - (b / a) * x1; 
    *y = x1; 
    return gcd; 
} 

int findGCD(int a, int b) {

    int x = 1, y = 1;
    return gcdExtended(a, b, &x, &y);
}

int modInverse(int e, int totient) {
    for(int d = 2; d < totient; d++) {
        if((e * d) % totient == 1) {
            return d;
        }
    }
    return -1;
}

// Function to compute base^expo mod m
int power(int base, int expo, int m) {
    long long res = 1;
    long long b = base % m;

    while(expo > 0) {
        if (expo & 1) {
            res = (res * b) % m;
        }
        b = (b * b) % m;
        expo = expo / 2;
    }
    return res;
}

// Encrypt message using public key (e, n)
int encrypt(int m, int e, int n) {
    return power(m, e, n);
}

// Decrypt message using private key (d, n)
int decrypt(int c, int d, int n) {
    return power(c, d, n);
}

int main() {
    
    int p = 1013;
    int q = 1019;

    int n = p * q;
    int totient = (p-1) * (q-1);
    
    int e = 3;  // e is the encrption exponent / public key
    int g = findGCD(e, totient);

    int d = modInverse(e, totient); // d is the decryption exponent / private key

    int M = 51010;
    int C = encrypt(M, e, n);
    int Mp = decrypt(C, d, n);

    printf("n             = %d\n", n);
    printf("totient       = %d\n", totient);
    printf("e             = %d\n", e);
    printf("d             = %d\n", d);
    printf("g             = %d\n", g);
     printf("e*d mod φ(n) = %d\n", e * d % totient);
    printf("M             = %d\n", M);
    printf("C=M^d%%n       = %lld\n", C);
    printf("M'            = %lld\n", Mp);
    printf("M == M'       ? %s\n", (Mp == M) ? "YES" : "NO");

    return 0;
}