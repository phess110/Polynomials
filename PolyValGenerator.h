#pragma once
#include <vector>
#include "Util.h"

class PolyValGenerator
{
private:
    std::vector<double> inputs;
    std::vector<double> outputs;
    std::vector<double> L;

public:
    PolyValGenerator(const std::vector<PtValPair> &);
    
    PolyValGenerator(const PolyValGenerator &) = delete;

    double Eval(double x) const;
};

