#pragma once
#include <vector>
#include "Util.h"

class PolyValGenerator
{
private:
    std::vector<double> m_inputs;
    std::vector<double> m_outputs;
    std::vector<double> m_L;

public:
    PolyValGenerator(const std::vector<PtValPair> &);
    
    PolyValGenerator(const PolyValGenerator &) = delete;

    double Eval(double x) const;
};

