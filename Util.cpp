#include <algorithm>
#include <complex>
#include "Util.h"

uint32_t pow2_round(uint32_t i) {
    return (i < 2) ? i :
        1 << static_cast<int> (std::ceil(std::log2(i)));
}

uint32_t bit_reverse(uint32_t v, uint32_t s) {
    uint32_t r = v & 1;
    s--;

    for (v >>= 1; v; v >>= 1)
    {
        r <<= 1;
        r |= v & 1;
        s--;
    }
    r <<= s;

    return r;
}

std::vector<cd> FFT(const std::vector<double> &a) {
    uint32_t N = a.size();
    uint32_t l = static_cast<int> (std::ceil(std::log2(N)));
    std::vector<cd> A(N);

    for (uint32_t k = 0; k < N; k++) {
        A[bit_reverse(k, l)] = a[k];
    }

    uint32_t m = 1;
    uint32_t n;
    cd wm, w, u, t;
    for (uint32_t s = 1; s <= l; s++) {
        n = m; // n = m / 2
        m <<= 1;
        wm = cd(std::cos(-TAU / m), std::sin(-TAU / m));
        for (uint32_t k = 0; k < N; k += m) {
            w = cd(1, 0);
            for (uint32_t j = 0; j <= n - 1; j++) {
                t = w * A[k + j + n];
                u = A[k + j];
                A[k + j] = u + t;
                A[k + j + n] = u - t;
                w *= wm;
            }
        }
    }

    return A;
}

std::vector<cd> InverseFFT(const std::vector<cd> &a) {
    uint32_t N = a.size();
    uint32_t l = static_cast<int> (std::ceil(std::log2(N)));
    std::vector<cd> A(N);

    for (uint32_t k = 0; k < N; k++) {
        A[bit_reverse(k, l)] = a[k];
    }

    uint32_t m = 1;
    uint32_t n;
    cd wm, w, u, t;
    for (uint32_t s = 1; s <= l; s++) {
        n = m; // n = m / 2
        m <<= 1;
        wm = cd(std::cos(TAU / m), std::sin(TAU / m));
        for (uint32_t k = 0; k < N; k += m) {
            w = cd(1, 0);
            for (uint32_t j = 0; j <= n - 1; j++) {
                t = w * A[k + j + n];
                u = A[k + j];
                A[k + j] = u + t;
                A[k + j + n] = u - t;
                w *= wm;
            }
        }
    }

    std::transform(A.begin(), A.end(), A.begin(), [N](cd x) { return x / static_cast<double>(N); });
    return A;
}