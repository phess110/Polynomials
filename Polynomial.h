#pragma once
#include <iostream> // cout
#include <type_traits>
#include <vector>

#include "Util.h"


/* Suggested TODOs
    - 
*/

class Polynomial
{
private:
    std::vector<double> m_coeffs;
    size_t m_degree;

    static std::vector<cd> PolyMultHelper(const Polynomial &, uint32_t);

    /*  Proper Usage:
        While input coefficients are allowed to be doubles, 
        it's best to only use this with integer coefficients. 

        For a polynomial with rational coefficients: f(x) = p0/q0 + ... + p(n)/q(n) x^n
        find the common denominator: d = lcm(q0,...,q(n)).
        Factor out 1/d to write f as a polynomial with integer coefficients: f(x) = 1/d * (a0 + ... + a(n)x^n).
        Perform desired operations on a0 + ... + a(n)x^n and then scale output by 1/d.
    */

    /* 
        Reverses the coefficients of the polynomial.
        Equivalent to computing x^n * p(1/x)
        param[in]: the polynomial to be reversed
    */
    static Polynomial ReversePolynomial(const Polynomial &); 

    /* Reverses coefficients of the polynomial in-place */
    void Reverse();

public:
    typedef std::pair<Polynomial, Polynomial> PolyPair;

    /* Constructor 
        Param [in]: vector A, the polynomial's coefficients: p(x) = sum a_n x^n.
    */
    Polynomial(const std::vector<double> &);

    /* Copy Constructor */
    Polynomial(const Polynomial &);

    /* Move Constructor */
    Polynomial(const Polynomial &&p) noexcept;

    /* Polynomial Multiplication via FFT
        param[in]: p, q
        param[in]: pow1, the power on p
        param[in]: pow2, the power on q
        return: the polynomial p(x)^pow1 * q(x)^pow2
    */
    static Polynomial PolyMult(const Polynomial &, const Polynomial &, 
                                uint8_t pow1 = 1, uint8_t pow2 = 1);

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
    static PolyPair PolyDiv(const Polynomial &, const Polynomial &);

    /*  Compute the Lagrange Coefficients 
        
        param[in]: A vector of point-value pairs {(x_i,y_i)}
        return: A vector L_i(x) = y_i * prod(j != i, 1/(x_i-x_j))

        P_i(x) = prod(j != i, (x-x_j)/(x_i-x_j))
        P_interp = sum(i=1,n+1, y_iP_i(x))

        For numerical stability, I don't optimize the calculation of the L_i's
    */
    static std::vector<double> PolyInterpolate(const std::vector<PtValPair> &);

    /* Returns the polynomial corresponding to the derivative of the input 
        param[in]: the polynomial to be differentiated
    */
    static Polynomial PolyDerivative(const Polynomial &);

    /*
        Returns the antiderivative of the given polynomial
        with 0 as the constant term.
    */
    static Polynomial PolyAntiDerivative(const Polynomial &);

    /* Polynomial operator overloads */
    Polynomial operator*(const double &d) const;
    Polynomial operator*(const Polynomial &p) const;
    Polynomial::PolyPair operator/(const Polynomial &q) const;
    Polynomial operator-(const Polynomial &q) const;
    Polynomial operator+(const Polynomial &q) const;
    double operator[](const size_t &i) const;

    /* Move assignment operator */
    Polynomial &operator=(Polynomial &&other) noexcept;

    /* Newton's Method for root finding. 
       param[in]: the initial guess for the root
       param[in]: (optional) tolerance. Default value: 1e-9
       param[in]: (optional) max_iters. Default value: 1000
    */
    double NewtonsMethod(double, 
                         double tolerance = 1e-6, 
                         uint32_t max_iters = 1e3) const;

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

    /*
        Computes the definite integral of the polynomial over the given interval.
    */
    double PolyIntegrate(double, double) const;

    /* Prints the polynomial to stdout */
    void PolyPrint() const;

};

