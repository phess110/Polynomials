#include <algorithm>
#include <chrono>         // std::chrono::milliseconds
#include <vector>

#include "PolyValGenerator.h"
#include "Polynomial.h"
#include "Util.h"

int main() {
    //auto t1 = std::chrono::high_resolution_clock::now();
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    //std::cout << duration << std::endl;
    //std::cout << std::this_thread::get_id() << std::endl;

    std::vector<PtValPair> v = { {1,2}, {2,3}, {5,7} };
    PolyValGenerator p = PolyValGenerator(v);

    std::cout << p.Eval(0.5) <<std::endl;
    std::cout << p.Eval(3) << std::endl;
    std::cout << p.Eval(4) << std::endl;
    std::cout << p.Eval(-10) << std::endl;

    //Polynomial::PolyPair r = Polynomial::PolyDiv(p, q);

    return 0;
}

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
