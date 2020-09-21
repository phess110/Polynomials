#pragma once
#include <iostream> // cout
#include <type_traits>
#include <vector>

#include "Util.h"


/* Suggested TODOs
    - 
*/



/*! 
    \class Polynomial

    \brief Represents a polynomial as a list of coefficients. 

    \remark Proper Usage: While input coefficients are allowed to be doubles, 
    it's best to only use this with integer coefficients.

    \remark For a polynomial with rational coefficients: \f$ f(x) = p_0/q_0 + ... + p_n/q_n x^n \f$,
    find the common denominator: \f$ d = \mathrm{lcm}(q_0,...,q_n) \f$.
    Factor out \f$ 1/d \f$ to write f as a polynomial with integer coefficients: \f$ f(x) = 1/d * (a_0 + ... + a_nx^n) \f$
    Perform desired operations on \f$ a_0 + ... + a_nx^n \f$ and then scale output by \f$ 1/d \f$.
*/

class Polynomial
{
private:
    std::vector<double> m_coeffs;
    size_t m_degree;

    static std::vector<cd> PolyMultHelper(const Polynomial &, uint32_t);

    /*! 
        \brief Equivalent to computing \f$ x^n * p(1/x) \f$
        \param [in] p the polynomial to be reversed
        \return the polynomial with the coefficients of p reversed
    */
    static Polynomial ReversePolynomial(const Polynomial &); 

    /*! Reverses coefficients of the polynomial in-place */
    void Reverse();

public:
    typedef std::pair<Polynomial, Polynomial> PolyPair;

    /*! Constructor 
        \param [in] A the vector of the polynomial's coefficients: \f$ p(x) = \sum a_n x^n \f$
    */
    Polynomial(const std::vector<double> &);

    /*! Copy Constructor */
    Polynomial(const Polynomial &);

    /*! Move Constructor */
    Polynomial(const Polynomial &&p) noexcept;

    /*! 
        \brief Polynomial Multiplication via FFT
        \param [in] p
        \param [in] q
        \param [in] pow1 the power of p. Optional. Default = 1
        \param [in] pow2 the power of q. Optional. Default = 1
        \return the polynomial \f$ p(x)^{pow1} * q(x)^{pow2} \f$
    */
    static Polynomial PolyMult(const Polynomial &, const Polynomial &, 
                                uint8_t pow1 = 1, uint8_t pow2 = 1);
    
    /*!
        \brief Compute inverse series of a polynomial
        \param [in] p the polynomial to be inverted
        \param [in] t a positive integer, the number of terms of its inverse to compute
        \return the power series of the inverse of the polynomial to desired number of terms 

        \warning the constant term p(0) MUST be non-zero.
        \throw std::invalid_argument Occurs when p(0) = 0
    */
    static Polynomial PolyInverse(const Polynomial &, uint32_t);

    /*! 
        \brief Polynomial division

        \param [in] f the dividend
        \param [in] g the divisor
        \return two polynomials, q(x) and r(x), such that \f$ f(x) = q(x)g(x) + r(x) \f$ 
     */
    static PolyPair PolyDiv(const Polynomial &, const Polynomial &);

    /*!  
        \brief Computes the Lagrange Coefficients 
        
        \param [in] points A vector of point-value pairs \f$ \{(x_i,y_i)\} \f$
        \return A vector \f$ L_i(x) = y_i  \prod_{j \neq i} \frac{1}{x_i-x_j} \f$ 

        \details Define \f$ P_i(x) = \prod_{j \neq i} \frac{x-x_j}{x_i-x_j} \f$. Then the interpolated
        polynomial is \f$ P_{interp} = \sum_{i=1}^{n+1} y_iP_i(x) \f$

        \note For numerical stability, I don't optimize the calculation of the L_i's
    */
    static std::vector<double> PolyInterpolate(const std::vector<PtValPair> &);

    /*!
        \param [in] p the polynomial to be differentiated
        \return the polynomial corresponding to the derivative of the input 
    */
    static Polynomial PolyDerivative(const Polynomial &);

    /*!
        \param [in] p the polynomial to integrate
        \return the antiderivative of p with 0 as the constant term.
    */
    static Polynomial PolyAntiDerivative(const Polynomial &);

    /* Polynomial operator overloads */

    /*! Polynomial-scalar multiplication */
    Polynomial operator*(const double &d) const;

    /*! Polynomial-polynomial multiplication */
    Polynomial operator*(const Polynomial &p) const;

    /*! Polynomial-polynomial division */
    Polynomial::PolyPair operator/(const Polynomial &q) const;

    /*! Polynomial-polynomial subtraction */
    Polynomial operator-(const Polynomial &q) const;

    /*! Polynomial-polynomial addition */
    Polynomial operator+(const Polynomial &q) const;

    /*! Polynomial coefficient indexing */
    double operator[](const size_t &i) const;

    /*! Move assignment operator */
    Polynomial &operator=(Polynomial &&other) noexcept;

    /*! 
       \brief Newton's Method for root finding. 
       \param [in] guess the initial guess for the root
       \param [in] tolerance Optional. Default value: 1e-6
       \param [in] max_iters Optional. Default value: 1000

       \return the calculated root
    */
    double NewtonsMethod(double, 
                         double tolerance = 1e-6, 
                         uint32_t max_iters = 1e3) const;

    /*!
        \brief Horner's method to evaluate polynomial at a point.

        \param [in] x the point to evaluate the polynomial at
        \return p(x)
    */
    double PolyEval(double) const;

    /*! 
        \brief Differentiates the polynomial in-place
    */
    void PolyDifferentiate();

    /*!
        \brief Computes the definite integral of the polynomial over the given interval.
        \param [in] s the starting point
        \param [in] e the end point
        \return the integral of the polynomial from s to e
    */
    double PolyIntegrate(double, double) const;

    /*! \brief Prints the polynomial to stdout */
    void PolyPrint() const;

};

