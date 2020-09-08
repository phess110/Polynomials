#include <vector>
#include "PolyValGenerator.h"
#include "Util.h"

PolyValGenerator::PolyValGenerator(const std::vector<PtValPair> &A) {
    // TODO initialize inputs, outputs, L 
}

double PolyValGenerator::Eval(double x) const {
    uint32_t N = inputs.size();
    double X = 1;
    double S = 0;
    // TODO Check if x is one of the original inputs

    // Otherwise
    for (uint32_t i = 0; i < N; i++) {
        X *= (x - inputs[i]);
        S += L[i] / (x - inputs[i]);
    }

    return X * S;
}
