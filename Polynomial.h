#pragma once
#include <iostream> // cout
#include <type_traits>
#include <vector>

#include "Util.h"

class Polynomial
{
private:
    std::vector<double> m_coeffs;
    size_t m_degree;

    static std::vector<cd> PolyMultHelper(const Polynomial &, uint32_t);

    /* TODO
        Reverses the coefficients of the polynomial.
        Equivalent to computing x^n * p(1/x)
    */
    // void ReversePolynomial(const Polynomial &); 
public:
    /* Constructor 
        Param [in]: vector A, the polynomial's coefficients: p(x) = sum a_n x^n.
    */
    Polynomial(const std::vector<double> &);

    /* Copy Constructor */
    Polynomial(const Polynomial &);

    /* Polynomial Multiplication via FFT
        param[in]: p, q
        return: the polynomial p(x) * q(x) 
    */
    static Polynomial PolyMult(const Polynomial &, const Polynomial &);

    /* Compute the power series of the inverse of the polynomial to desired number of terms 
        param[in]: polynomial p to be inverted
        param[in]: positive integer t, the number of terms of its inverse to compute
        
        Note: the constant term p(0) MUST be non-zero. Otherwise, no inverse exists and the function aborts.
    */
    static Polynomial PolyInverse(const Polynomial &, uint32_t);

    /* Polynomial Division

        param[in]: f, g 
        return: two polynomials, q(x) and r(x), such that f(x) = q(x)g(x) + r(x) 
     */
    static std::pair<Polynomial, Polynomial> PolyDiv(const Polynomial &, const Polynomial &);

    /*  Compute the Lagrange Coefficients 
        
        param[in]: A vector of point-value pairs {(x_i,y_i)}
        return: A vector L_i(x) = y_i * prod(j != i, 1/(x_i-x_j))

        P_i(x) = prod(j != i, (x-x_j)/(x_i-x_j))
        P_interp = sum(i=1,n+1, y_iP_i(x))

        For numerical stability, I don't optimize the calculation of the L_i's
    */
    static std::vector<double> PolyInterpolate(const std::vector<PtValPair> &);

    /* Returns the polynomial corresponding to the derivative of the input */
    static Polynomial PolyDerivative(const Polynomial &);

    /* Newton's Method for root finding. 
       param[in]: the initial guess for the root
       param[in]: (optional) tolerance. Default value: 1e-9
       param[in]: (optional) max_iters. Default value: 1000
    */
    double NewtonsMethod(double, 
                         double tolerance = 1e-9, 
                         uint32_t max_iters = 1e3);

    /*
        Horner's method to evaluate polynomial at x.

        param[in]: double x 
        return: p(x)
    */
    double PolyEval(double) const;

    /* 
        Differentiates the polynomial
    */
    void PolyDifferentiate();

    /* Prints the polynomial to stdout */
    void Print() const;

};

