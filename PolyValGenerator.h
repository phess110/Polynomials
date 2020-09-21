#pragma once
#include <vector>
#include "Util.h"

/*!
    \class PolyValGenerator
    \brief The purpose of this class is to provide an efficient way of evaluating 
    an interpolated polynomial without actually having to compute the coefficients
    of that polynomial.
*/
class PolyValGenerator
{
private:
    std::vector<double> m_inputs;
    std::vector<double> m_outputs;
    std::vector<double> m_L;

public:
    /*! Constructs a generator using Polynomial::PolyInterpolate to compute the coefficients */
    PolyValGenerator(const std::vector<PtValPair> &);
    
    /* Remove copy constructor */
    PolyValGenerator(const PolyValGenerator &) = delete;

    /*! Move Constructor */
    PolyValGenerator(const PolyValGenerator &&) noexcept;

    /*! Evaluates the generator at the point x */
    double Eval(double x) const;
};