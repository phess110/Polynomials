#include <algorithm>
#include <chrono>         // std::chrono::milliseconds
#include <vector>

#include "PolyValGenerator.h"
#include "Polynomial.h"
#include "Util.h"


PolyValGenerator::PolyValGenerator(const std::vector<PtValPair> &A) {
    uint32_t N = A.size();
    for (uint32_t i = 0; i < N; i++) {
        PtValPair p = A[i];
        m_inputs.push_back(p.x);
        m_outputs.push_back(p.y);
    }
    m_L = Polynomial::PolyInterpolate(A);
}

PolyValGenerator::PolyValGenerator(const PolyValGenerator &&p) noexcept :
    m_inputs(std::move(p.m_inputs)),
    m_outputs(std::move(p.m_outputs)),
    m_L(std::move(p.m_L))
{
    std::cout << "Move construct\n";
}


double PolyValGenerator::Eval(double x) const {
    uint32_t N = m_inputs.size();
    double X = 1;
    double S = 0;

    auto it = std::find(m_inputs.begin(), m_inputs.end(), x);
    // Check if x is one of the original inputs
    if (it != m_inputs.end()) {
        int d = std::distance(m_inputs.begin(), it);
        return m_outputs[d];
    }

    // Otherwise
    for (uint32_t i = 0; i < N; i++) {
        X *= (x - m_inputs[i]);
        S += m_L[i] / (x - m_inputs[i]);
    }

    return X * S;
}
