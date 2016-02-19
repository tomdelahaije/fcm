/** Utilities to compute Spherical-Harmonics (SH)-related functions. Take a look to the
 corresponding implementation (.cxx) file for details on the meaning and usage of each
 of the functions in this file.
 
 IMPORTANT NOTE: Part of the implementations of the functions in this headers file is
 in the file sh2hot.cxx. Hence, it is mandatory to include sh2hot.cxx in any project
 including "sphericalHarmonics.h"
 */

#ifndef _sphericalHarmonics_h_
#define _sphericalHarmonics_h_

#include <vnl/vnl_matrix.h> // Use the powerful vnl routines for matrix manipulation
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_convolve.h>
#include <vector>

namespace shmaths
{
    
#ifndef PI
#define PI 3.14159265358979323846
#endif
    
#ifndef GAMMAEULER
#define GAMMAEULER 0.5772156649015328606065
#endif
    
#ifndef ROOTSOLVERTOL
#define ROOTSOLVERTOL 1e-6
#endif
    
#ifndef ROOTSOLVERTOLNULLCOEFFS
#define ROOTSOLVERTOLNULLCOEFFS 1e-12
#endif
    
#ifndef NULLPOLYNOMIALTHRESHOLD
#define NULLPOLYNOMIALTHRESHOLD 1e-4
#endif
    
#ifndef CANCELSCALEPOLYNOMIALSUBTRACT
#define CANCELSCALEPOLYNOMIALSUBTRACT 1e3
#endif
    
    
    typedef vnl_matrix<double>                  SHMatrixType;
    typedef vnl_vector<double>                  PolynomialType;
    typedef vnl_vector< vcl_complex<double> >   RootsType;
    typedef vnl_matrix<unsigned int>            DegreesType;
    
    unsigned long cumprod( const unsigned int start, const unsigned int end );
    
    void computeAssociatedLegendrePolynomials( const double x, const unsigned int L, double* buffer );
    
    void computeLegendrePolynomials( const double x, const unsigned int L, double* buffer );
    
    double* allocateBufferForAssociatedLegendrePolynomials( const unsigned int L );
    
    double* allocateBufferForLegendrePolynomials( const unsigned int L );
    
    void computeEvenAssociatedLegendrePolynomials( const double x, const unsigned int L, double* buffer );
    
    void computeEvenLegendrePolynomials( const double x, const unsigned int L, double* buffer );
    
    double* allocateBufferForEvenAssociatedLegendrePolynomials( const unsigned int L );
    
    double* allocateBufferForEvenLegendrePolynomials( const unsigned int L );
    
    unsigned int getNumberOfLegendrePolynomials( const unsigned int L );
    
    unsigned int getNumberOfAssociatedLegendrePolynomials( const unsigned int L );
    
    unsigned int getNumberOfEvenLegendrePolynomials( const unsigned int L );
    
    unsigned int getNumberOfEvenAssociatedLegendrePolynomials( const unsigned int L );
    
    double computeP_l( const double x, const unsigned int l );
    
    double computeP_l_m( const double x, const unsigned int l, const int m );
    
    double computeP_l( const unsigned int l, const double* buffer );
    
    double computeP_l_m( const unsigned int l, const int, const double* buffer );
    
    double computeY_l_m( const double phi, const double theta, const unsigned int l, const int m );
    
    void computeSHMatrix( const unsigned int N, const double* theta, const double* phi, const unsigned int L, SHMatrixType& sh  );
    
    void computeSHMatrixSymmetric( const unsigned int N, const double* theta, const double* phi, const unsigned int L, SHMatrixType& sh  );
    
    void computeSHEigMatrix( const unsigned int L, SHMatrixType& eig );
    
    void computeSHEigMatrixSymmetric( const unsigned int L, SHMatrixType& eig );
    
    void computeSHFRTMatrix( const unsigned int L, SHMatrixType& frt );
    
    void computeSHFRTMatrixSymmetric( const unsigned int L, SHMatrixType& frt );
    
    double expint( const double x );
    
    double Ein( const double x );
    
    void computeShericalCoordsFromCartesian( const double x, const double y, const double z, double& theta, double& phi );
    
    void computeShericalCoordsFromCartesian( const double* x, const double* y, const double* z, double* theta, double* phi, const unsigned int N );
    
    /** THE FOLLOWING CODE IS USED TO TRANSLATE BETWEEN SPHERICAL HARMONICS AND HIGHER ORDER TENSORS.
     These functions are implemented in the file sh2hot.cxx */ 
    
    void computeSetOfHOTPowers( const unsigned int order, unsigned int* nx, unsigned int* ny, unsigned int* nz );
    
    void computeMultiplicityOfHOTComponents( const unsigned int order, unsigned long* mu, const unsigned int* nx, const unsigned int* ny, const unsigned int* nz );
    
    bool generateSH2HOTMatrix( const unsigned int order, SHMatrixType& M );
    
    void evaluateHOTFromCoefficients( const unsigned int N, const double* x, const double* y, const double* z, double* val, const unsigned int order,
                                     const unsigned int* nx, const unsigned int* ny, const unsigned int* nz, const unsigned long* mu, const SHMatrixType T );
    
    SHMatrixType computeReducedTensor( const SHMatrixType original, const unsigned int order, const unsigned int coordinate );
    
    
    void computeLocalExtrema( const unsigned int order, const SHMatrixType hot, std::vector< vnl_vector<double> >& output );
    
    void computeSolutionForm100( const unsigned int order, const SHMatrixType hotx, const SHMatrixType hoty, const SHMatrixType hotz,
                                std::vector< vnl_vector<double> >& output );
    
    void computeSolutionFormt10( const unsigned int order, const SHMatrixType hotx, const SHMatrixType hoty, const SHMatrixType hotz, 
                                std::vector< vnl_vector<double> >& output, const unsigned long* mu );
    
    void computeSolutionFormtu1( const unsigned int order, const SHMatrixType hotx, const SHMatrixType hoty, const SHMatrixType hotz,
                                std::vector< vnl_vector<double> >& output, const unsigned long* mu );
    
    int polynomialRoots( const PolynomialType& polynomial, RootsType& roots );
    
    int polynomialSystem( const PolynomialType& P, const PolynomialType& Q, const DegreesType& degP,
                         const DegreesType& degQ, RootsType& rootsX, RootsType& rootsY );
    
    bool polynomialIsNull( const PolynomialType& P );
    
    double normalizePolynomial( std::vector<PolynomialType>& P, unsigned int degree );
    
    double normalizePolynomial( PolynomialType& polynomial );
    
    double evaluatePolynomial( const PolynomialType& polynomial, const double x );
    
    PolynomialType evaluatePolynomial( const std::vector<PolynomialType>& polynomial, const double y );
    
    PolynomialType addPolynomials( const PolynomialType& pol1, const PolynomialType& pol2 );
    
    PolynomialType polynomialGCD( const PolynomialType& pol1, const PolynomialType& pol2 );
    
    PolynomialType polynomialDivision( const PolynomialType& pol1, const PolynomialType& pol2, PolynomialType& rem );
    
    
} // End namespace shmaths

#endif //_sphericalHarmonics_h_
