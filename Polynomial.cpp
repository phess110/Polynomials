#include <vector>

#include "Polynomial.h"

#include <iostream>       // std::cout
#include <future>         // std::async, std::future
#include <chrono>         // std::chrono::milliseconds

// async test
bool is_prime(int x) {
    std::cout << std::this_thread::get_id() << std::endl;
    for (int i = 2; i < x; ++i) if (x % i == 0) return false;
    return true;
}

// test async stuff
void test() {  
    std::future<bool> fut = std::async(std::launch::async, is_prime, 1000000241);
    std::future<bool> fut2 = std::async(std::launch::async, is_prime, 1000000207);
    std::cout << "checking...\n";

    fut.wait();
    fut2.wait();

    std::cout << fut.get() << std::endl;
    std::cout << fut2.get() << std::endl;
}

int main() {
    //std::cout << std::this_thread::get_id() << std::endl;
    //auto t1 = std::chrono::high_resolution_clock::now();
    //test();
    //auto t2 = std::chrono::high_resolution_clock::now();

    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    //std::cout << duration << std::endl;
    //std::cout << std::this_thread::get_id() << std::endl;

    Polynomial p = Polynomial({ 0,0,0,1 });
    Polynomial q = Polynomial({ 0,0,0,1 });
    Polynomial r = Polynomial::PolyMult(p, q);
    return 0;
}

Polynomial::Polynomial(const std::vector<double> &A) : m_degree(A.size() - 1) { // TODO the vector should be nonempty!
    m_coeffs = std::vector<double>(A);
}

Polynomial::Polynomial(const Polynomial &p) : m_degree(p.m_degree) {
    m_coeffs = std::vector<double>(p.m_coeffs);
}

std::vector<cd> Polynomial::PolyMultHelper(const Polynomial &p, uint32_t N) {
    std::vector<double> p_coeffs(p.m_coeffs);
    
    while (p_coeffs.size() < N) {
        p_coeffs.push_back(0);
    }

    return FFT(p_coeffs);
}

Polynomial Polynomial::PolyMult(const Polynomial &p, const Polynomial &q) {
    uint32_t N = pow2_round(p.m_degree + q.m_degree + 1);

    // Could be parallelized
    std::vector<cd> pFFT = PolyMultHelper(p, N);
    std::vector<cd> qFFT = PolyMultHelper(q, N);
    
    // Take the Hadamard product
    for (uint32_t i = 0; i < N; i++) {
        pFFT[i] *= qFFT[i];
    }

    // Use IFFT to recover product coeffs
    std::vector<cd> out = InverseFFT(pFFT);
    std::vector<double> out_real(N);
    std::transform(out.begin(), out.end(), out_real.begin(), [](cd x) { return std::real(x); });
    return Polynomial(out_real);
}

Polynomial Polynomial::PolyInverse(const Polynomial &p, uint32_t t) {
    uint32_t m = 1;

    // TODO check leading coeff is nonzero

    Polynomial inv = Polynomial({ 1 / p.m_coeffs[0] });
    while (m < t) {
        m <<= 1;
        // TODO 
        // inv = 2 * inv - A * inv^2
    }

    return inv;
}

std::pair<Polynomial, Polynomial> Polynomial::PolyDiv(const Polynomial &f, const Polynomial &g) {

    return std::pair<Polynomial, Polynomial>(f, g); // TODO 
}

std::vector<double> Polynomial::PolyInterpolate(const std::vector<PtValPair> &points) {
    uint32_t N = points.size();
    std::vector<double> L(N);

    double l;
    uint32_t i, j;
    for (i = 0; i < N; i++) {
        l = points[i].y;
        for (j = 0; j < N; j++) {
            if (j != i) {
                l /= (points[i].x - points[j].x);
            }
        }
        L[i] = l;
    }

    return L;
}

Polynomial Polynomial::PolyDerivative(const Polynomial &p) {
    std::vector<double> deriv(p.m_degree);
    for (uint32_t i = 0; i < p.m_degree; i++) {
        deriv[i] = (static_cast<double>(i) + 1) * p.m_coeffs[i + 1];
    }

    return Polynomial(deriv); 
}

double Polynomial::PolyEval(double x) const {
    double p = m_coeffs[m_degree] * x;
    for (uint32_t i = m_degree - 1; i > 0; i--) {
        p += m_coeffs[i];
        p *= x;
    }
    return m_coeffs[0] + p;
}

void Polynomial::PolyDifferentiate() {
    for (uint32_t i = 0; i < m_degree; i++) {
        m_coeffs[i] = (static_cast<double>(i) + 1) * m_coeffs[i + 1];
    }
    m_coeffs[m_degree] = 0;
}

double Polynomial::NewtonsMethod(double guess, double tolerance, uint32_t max_iters) {
    Polynomial deriv = PolyDerivative(*this);

    double x0 = guess, x1;

    uint32_t i = 0;
    while (i < max_iters)
    {
        x1 = x0 - PolyEval(x0) / deriv.PolyEval(x0);

        if (std::abs(x1 - x0) < tolerance) {
            break;
        }

        i++;
        x0 = x1;
    }

    return x0;
}


void Polynomial::Print() const {
    for (uint32_t i = 0; i < m_degree; i++) {
        std::cout << m_coeffs[i] << "x^" << i << " + ";
    }
    std::cout << m_coeffs[m_degree] << "x^" << m_degree;
}
