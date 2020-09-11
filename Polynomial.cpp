
#include <future>         // std::async, std::future
#include <iostream>       // std::cout
#include <vector>

#include "Polynomial.h"
#include "Util.h"

#include <chrono>         // std::chrono::milliseconds

int main() {
    //auto t1 = std::chrono::high_resolution_clock::now();
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    //std::cout << duration << std::endl;
    //std::cout << std::this_thread::get_id() << std::endl;

    Polynomial p = Polynomial({1,1,1,1,-1});
    //Polynomial q = Polynomial({1,2,3,4,5,6,7});
    p.Print();
    //std::cout << "*\n";
    //q.Print();
    //std::cout << "=\n";
    Polynomial r = Polynomial::PolyInverse(p, 8);
    r.Print();
    return 0;
}

Polynomial::Polynomial(const std::vector<double> &A) : m_degree(A.size() - 1) { // TODO the vector should be nonempty!
    m_coeffs = std::vector<double>(A);
}

Polynomial::Polynomial(const Polynomial &p) : m_degree(p.m_degree) {
    m_coeffs = std::vector<double>(p.m_coeffs);
    std::cout << "Copied\n";
}

std::vector<cd> Polynomial::PolyMultHelper(const Polynomial &p, uint32_t N) {
    std::vector<double> p_coeffs(p.m_coeffs);
    
    while (p_coeffs.size() < N) {
        p_coeffs.push_back(0);
    }

    return FFT(p_coeffs);
}

Polynomial Polynomial::PolyMult(const Polynomial &p, const Polynomial &q,
                                uint8_t pow1, uint8_t pow2) {
    uint32_t num_coeffs = pow1 * p.m_degree + pow2 * q.m_degree + 1;
    uint32_t N = pow2_round(num_coeffs);

    // Parallelize FFT computation on p and q
    std::future<std::vector<cd>> f1 = std::async(PolyMultHelper, p, N);
    std::future<std::vector<cd>> f2 = std::async(PolyMultHelper, q, N);
    f1.wait();
    f2.wait();
    std::vector<cd> pFFT = f1.get();
    std::vector<cd> qFFT = f2.get();
    
    // Compute r = p^pow1 * q^pow2 on primitive N-th roots of unity
    for (uint32_t i = 0; i < N; i++) {
        pFFT[i] = std::pow(pFFT[i], pow1) * std::pow(qFFT[i], pow2);
    }

    // Use IFFT to recover product coeffs of r
    std::vector<cd> out = InverseFFT(pFFT);
    std::vector<double> out_real(num_coeffs);
    /*
        Since we started with a real polynomial, the IFFT should have no imaginary part.
        Due to lack of precision, the imaginary part may be nonzero, but we can safely ignore it.

        We need a reliable way of hiding the imprecision in the real part. 
        Currently, I just add and then subtract 10 to truncate any really small errors. 
    */
    std::transform( out.begin(), 
                    out.begin() + num_coeffs, 
                    out_real.begin(), 
                    [](cd x) { return roundError(std::real(x)); });
    return Polynomial(out_real);
}

Polynomial Polynomial::PolyInverse(const Polynomial &p, uint32_t t) {
    uint32_t m = 1;

    // TODO throw exception???
    if (p.m_coeffs[0] == 0) {
        // error

        return Polynomial({ 0 });
    }

    Polynomial inv = Polynomial({ 1 / p.m_coeffs[0] });
    while (m < t) {
        m <<= 1;
        // inv = 2 * inv - A * inv^2
        inv = (inv * 2) - PolyMult(p, inv, 1, 2);
        inv.m_coeffs.resize(m);
        inv.m_degree = m - 1;
    }

    return inv;
}

/*
    TODO
    Update all functions to ignore terms with too high degree
*/

Polynomial Polynomial::operator*(const double& d) {
    std::vector<double> result(m_degree + 1);

    std::transform( m_coeffs.begin(),
                    m_coeffs.end(),
                    result.begin(),
                    [d](double x) { return d * x; });
    return Polynomial(result);
}

Polynomial Polynomial::operator*(const Polynomial &p) {
    return PolyMult(*this, p);
}

Polynomial::PolyPair Polynomial::operator/(const Polynomial &q) {
    return PolyDiv(*this, q);
}

Polynomial Polynomial::operator-(const Polynomial &q) {
    uint32_t degree = std::max(m_degree, q.m_degree);

    std::vector<double> result(degree + 1);
    double c, d;
    for (uint32_t i = 0; i < degree + 1; i++) {
        c = (m_degree < i) ? 0 : m_coeffs[i];
        d = (q.m_degree < i) ? 0 : q.m_coeffs[i];
        result[i] = c - d;
    }

    return Polynomial(result);
}

Polynomial Polynomial::operator+(const Polynomial &q) {
    uint32_t degree = std::max(m_degree, q.m_degree);

    std::vector<double> result(degree + 1);
    double c, d;
    for (uint32_t i = 0; i < degree + 1; i++) {
        c = (m_degree < i) ? 0 : m_coeffs[i];
        d = (q.m_degree < i) ? 0 : q.m_coeffs[i];
        result[i] = c + d;
    }

    return Polynomial(result);
}

Polynomial::Polynomial(const Polynomial &&p) noexcept : 
    m_degree(std::move(p.m_degree)),
    m_coeffs(std::move(p.m_coeffs)) { std::cout << "Move construct\n"; }

Polynomial & Polynomial::operator=(Polynomial &&other) noexcept {
    if (this != &other) {
        m_degree = std::move(other.m_degree);
        m_coeffs = std::move(other.m_coeffs);
    }
    std::cout << "Move assign\n";
    return *this;
}

// TODO implement operator overloads
/*
    1. Scalar multiplication
    2. Polynomial addition
    3. Polynomial subtraction
    4. Polynomial multiplication
    5. Polynomial division
    6. Move constructor
    7. Move assignment op
*/

Polynomial::PolyPair Polynomial::PolyDiv(const Polynomial &f, const Polynomial &g) {
    // TODO
/*
    uint32_t N = f.m_degree - g.m_degree + 1;
    Polynomial fR = PolyReverse(f);
    Polynomial gR = PolyReverse(g);

    Polynomial qR = PolyMult(fR, PolyInverse(gR, ))

    Polynomial q = reverse(qR)
    Polynomial r = f - q * g;
    return std::pair<Polynomial, Polynomial>(q, r);
*/

    return std::pair<Polynomial, Polynomial>(f, g);
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
    // TODO update m_coeff size  -> remove last entry
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
    std::cout << m_coeffs[m_degree] << "x^" << m_degree << std::endl;
}
