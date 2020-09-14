
#include <future>         // std::async, std::future
#include <iostream>       // std::cout
#include <vector>

#include "Polynomial.h"
#include "Util.h"

Polynomial::Polynomial(const std::vector<double> &A) { 
    if (A.empty()) {
        m_degree = 0;
        m_coeffs = { 0 };
    }
    else {
        m_degree = A.size() - 1;
        m_coeffs = std::vector<double>(A);
    }
}

Polynomial::Polynomial(const Polynomial &p) : m_degree(p.m_degree) {
    m_coeffs = std::vector<double>(p.m_coeffs);
    std::cout << "Copied\n";
}

Polynomial::Polynomial(const Polynomial &&p) noexcept :
    m_degree(std::move(p.m_degree)),
    m_coeffs(std::move(p.m_coeffs)) {
    std::cout << "Move construct\n";
}

Polynomial &Polynomial::operator=(Polynomial &&other) noexcept {
    if (this != &other) {
        m_degree = std::move(other.m_degree);
        m_coeffs = std::move(other.m_coeffs);
    }
    std::cout << "Move assign\n";
    return *this;
}

Polynomial Polynomial::ReversePolynomial(const Polynomial &p) {
    uint32_t N = p.m_degree + 1;
    std::vector<double> rev(N);
    uint32_t endZeros = 0;
    double d;
    for (uint32_t i = 0; i < N; i++) {
        d = p[N - i - 1];
        rev[i] = d;
        if (d == 0) {
            endZeros++;
        }
        else {
            endZeros = 0;
        }
    }

    rev.resize(N - endZeros);
    return Polynomial(rev);
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
    std::future<std::vector<cd>> f1 = std::async(PolyMultHelper, std::cref(p), N);
    std::future<std::vector<cd>> f2 = std::async(PolyMultHelper, std::cref(q), N);
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
        Due to lack of precision, the imaginary part may be nonzero, so we can safely ignore it.

        Use some rounding to hide the imprecision in the real part. 
    */
    std::transform( out.begin(), 
                    out.begin() + num_coeffs, 
                    out_real.begin(), 
                    [](cd x) { return roundError(std::real(x)); });
    return Polynomial(out_real);
}

Polynomial Polynomial::PolyInverse(const Polynomial &p, uint32_t t) {
    if (p[0] == 0) {
        throw "Inverse does not exist";
    }

    uint32_t m = 1;
    Polynomial inv = Polynomial({ 1 / p[0] });
    while (m < t) {
        m <<= 1;
        // inv = 2 * inv - A * inv^2
        inv = (inv * 2) - PolyMult(p, inv, 1, 2);
        inv.m_coeffs.resize(m);
        inv.m_degree = m - 1;
    }
    return inv;
}

Polynomial Polynomial::operator*(const double& d) const {
    if (d == 0) {
        return Polynomial({ 0 });
    }
    else {
        std::vector<double> result(m_degree + 1);

        std::transform(m_coeffs.begin(),
            m_coeffs.end(),
            result.begin(),
            [d](double x) { return d * x; });
        return Polynomial(result);
    }
}

Polynomial Polynomial::operator*(const Polynomial &p) const {
    return PolyMult(*this, p);
}

Polynomial::PolyPair Polynomial::operator/(const Polynomial &q) const {
    return PolyDiv(*this, q);
}

Polynomial Polynomial::operator-(const Polynomial &q) const {
    uint32_t num_coeffs = std::max(m_degree, q.m_degree) + 1;

    std::vector<double> result(num_coeffs);
    double c, d;
    uint8_t trailingZeros = 0;
    for (uint32_t i = 0; i < num_coeffs; i++) {
        c = (m_degree < i) ? 0 : m_coeffs[i];
        d = (q.m_degree < i) ? 0 : q[i];
        result[i] = c - d;
        if (result[i] == 0) {
            trailingZeros++;
        }
        else {
            trailingZeros = 0;
        }
    }

    result.resize(num_coeffs - trailingZeros);
    return Polynomial(result);
}

Polynomial Polynomial::operator+(const Polynomial &q) const {
    uint32_t num_coeffs = std::max(m_degree, q.m_degree) + 1;

    std::vector<double> result(num_coeffs);
    double c, d;
    uint8_t trailingZeros = 0;
    for (uint32_t i = 0; i < num_coeffs; i++) {
        c = (m_degree < i) ? 0 : m_coeffs[i];
        d = (q.m_degree < i) ? 0 : q[i];
        result[i] = c + d;
        if (result[i] == 0) {
            trailingZeros++;
        }
        else {
            trailingZeros = 0;
        }
    }

    result.resize(num_coeffs - trailingZeros);
    return Polynomial(result);
}

double Polynomial::operator[](const size_t &i) const {
    return m_coeffs[i];
}

Polynomial::PolyPair Polynomial::PolyDiv(const Polynomial &f, const Polynomial &g) {
    uint32_t N = f.m_degree - g.m_degree + 1;
    Polynomial fR = ReversePolynomial(f);
    Polynomial gR = ReversePolynomial(g);

    Polynomial qR = PolyMult(fR, PolyInverse(gR, N));
    qR.m_coeffs.resize(N);
    qR.m_degree = N - 1;
    qR.Reverse();

    Polynomial r = f - (qR * g);

    return PolyPair(qR, r);
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
    if (p.m_degree == 0) {
        return Polynomial({ 0 });
    }
    else {
        std::vector<double> deriv(p.m_degree);
        for (uint32_t i = 0; i < p.m_degree; i++) {
            deriv[i] = (static_cast<double>(i) + 1) * p[i + 1];
        }

        return Polynomial(deriv);
    }
}

double Polynomial::PolyEval(double x) const {
    if (m_degree == 0) {
        return m_coeffs[0];
    }
    else {
        double p = m_coeffs[m_degree] * x;
        for (uint32_t i = m_degree - 1; i > 0; i--) {
            p += m_coeffs[i];
            p *= x;
        }
        return m_coeffs[0] + p;
    }
}

void Polynomial::PolyDifferentiate() {
    for (uint32_t i = 0; i < m_degree; i++) {
        m_coeffs[i] = (static_cast<double>(i) + 1) * m_coeffs[i + 1];
    }
    m_coeffs[m_degree] = 0;
    if (m_degree > 0) {
        m_degree--;
        m_coeffs.pop_back();
    }
}

double Polynomial::NewtonsMethod(double guess, 
                                double tolerance, 
                                uint32_t max_iters) const {
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
    std::cout.precision(15);
    std::cout << "After " << i << " iterations: " << x0 << std::endl;
    return x0;
}

void Polynomial::Reverse() {
    std::reverse(std::begin(m_coeffs), std::end(m_coeffs));
    int i;
    for (i = m_coeffs.size() - 1; i >= 1 && m_coeffs[i] == 0; i--);
    m_coeffs.resize(i + 1);
    m_degree = i;
}


void Polynomial::Print() const {
    for (uint32_t i = 0; i < m_degree; i++) {
        std::cout << m_coeffs[i] << "x^" << i << " + ";
    }
    std::cout << m_coeffs[m_degree] << "x^" << m_degree << std::endl;
}
