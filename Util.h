#pragma once
#include <complex>
#include <vector>

/*
    Provides some useful functionality that's not inherently tied to the Polynomial class.
*/

/*
    Suggested TODO
    Use FFT/Convolution to implement string pattern matching (with/without wildcards)
*/

typedef std::complex<double> cd;
const double PI = 3.1415926535897932384626;
const double TAU = 6.283185307179586476925;
const double epsilon = 1e-6;

typedef struct PtValPair {
    double x;
    double y;
} PtValPair;

/* Return the ceil(x/2) of integer x */
#define ceildiv2(x) ((x >> 1) + (x & 1));

/* Rounds i up to nearest power of 2 */
uint32_t pow2_round(uint32_t i);

/* If x is within epsilon of an integer, round x */
double roundError(double x);

/* 
    Reverse the bits of v. 
    v will be treated an an s-bit number.
*/
uint32_t bit_reverse(uint32_t v, uint32_t s);

/* 
    Iterative Fast Fourier Transform 
    Note: The size(a) should be a power of 2.
*/
std::vector<cd> FFT(const std::vector<double> &a);

/* 
    Iterative Inverse Fast Fourier Transform 
    Note: The size(a) should be a power of 2.
 */
std::vector<cd> InverseFFT(const std::vector<cd> &a);