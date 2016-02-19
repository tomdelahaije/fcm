/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkComputeLocalCostFromSH.txx,v $
 Language:  C++
 Date:      $Date: 2006/01/11 19:43:31 $
 Version:   $Revision: 1.21 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef _itkComputeLocalCostFromSH_txx
#define _itkComputeLocalCostFromSH_txx

#include "itkComputeLocalCostFromSH.h"
#include "itkImageRegionIterator.h"
#include "math.h"
// This is used only in the MEANTIME case; this header files
// contains a function to compute the scaling factors to sharpen
// the ODF as a function of the diffusivity of the corresponding
// voxel:
#include "GenerateDeconvolutionFactors.h"


//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
/*
#ifndef CODE_FOR_DEBUG_EXPERIMENTAL
#define CODE_FOR_DEBUG_EXPERIMENTAL
#endif
*/
// DEBUG
#ifdef CODE_FOR_DEBUG_EXPERIMENTAL
#include "itkImageFileWriter.h"
#endif
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

namespace itk
{
    template< class TInputImage, class TOutputImage >
    ComputeLocalCostFromSH< TInputImage, TOutputImage >
    ::ComputeLocalCostFromSH()
    {
        if( TInputImage::ImageDimension != 3 )
            itkExceptionMacro( << "This filter is only supported for input dimension 3" );
        m_NeighboringDirections.resize(0);
        m_LocalCostType = ATTODF;
        m_SquaredCosines.SetSize(0,0);
        
        m_bD.clear();
        m_LUT.clear();
        
        m_L = 6;
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
#ifdef CODE_FOR_DEBUG_EXPERIMENTAL
        // DEBUG
        this->SetNumberOfThreads(1);
#endif
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
    }
    
    template< class TInputImage, class TOutputImage >
    void ComputeLocalCostFromSH< TInputImage, TOutputImage >
    ::GenerateOutputInformation( void )
    {
        Superclass::GenerateOutputInformation();
        this->GetOutput()->SetVectorLength( m_NeighboringDirections.size() );
    }
    
    template< class TInputImage, class TOutputImage >
    void ComputeLocalCostFromSH< TInputImage, TOutputImage >
    ::BeforeThreadedGenerateData( void )
    {
        // Call superclass' implementation for this method:
        Superclass::BeforeThreadedGenerateData();
        // Compute a couple of parameters of interest and check:
        unsigned int N = m_NeighboringDirections.size();
        if( N==0 )
            itkExceptionMacro( << "You must set the directions to compute the costs before executing the filter" );
        unsigned int M = this->GetInput()->GetVectorLength();
        // From M, we need to compute the order of the SH:
        switch(M){
            case 1:
                m_L = 0;
                break;
            case 6:
                m_L = 2;
                break;
            case 15:
                m_L = 4;
                break;
            case 28:
                m_L = 6;
                break;
            case 45:
                m_L = 8;
                break;
            default:
                itkExceptionMacro( << "Only degrees 0, 2, 4, 6, and 8 are supported. (There are " << M << " coefficients)" );
                break;
        }
        
        // The directions to compute the costs are given in cartesian coordinates. We
        // need to use spherical coordinates in order to compute the matrix to convert
        // from SH coefficients to actual costs. This is done by the shmaths code:
        double* gx = new double[N]; // Old-style c vector
        double* gy = new double[N]; // Old-style c vector
        double* gz = new double[N]; // Old-style c vector
        // Retrieve the cartesian coordinates in standard c vectors:
        for( unsigned int k=0; k<N; ++k ){
            if( m_NeighboringDirections[k].GetNorm()>1e-6 )
                m_NeighboringDirections[k].Normalize();
            gx[k] = m_NeighboringDirections[k][0];
            gy[k] = m_NeighboringDirections[k][1];
            gz[k] = m_NeighboringDirections[k][2];
        }
        // Now, compute the spherical coordinates corresponding to these gradient directions:
        double* theta = new double[N];
        double* phi   = new double[N];
        shmaths::computeShericalCoordsFromCartesian( gx, gy, gz, theta, phi, N );
        // Cartesian coordinates are no longer necessary:
        delete[] gx; delete[] gy; delete[] gz;
        // From the spherical coordinates, we can compute all required matrixes (standard convention):
        shmaths::computeSHMatrixSymmetric( N, theta, phi, m_L, m_B  ); // B has size N x M
        shmaths::computeSHFRTMatrixSymmetric( m_L, m_FRT );            // FRT matrix has size M x M
        // Use the appropriate matrices depending on how we compute the local costs:
        MatrixType eigenMatrix;
        switch( m_LocalCostType ){
            case ATTODF:
                m_FRT = m_B*m_FRT;
                break;
            case MEANTIME:
                /** In this case, we have to weight the costs in this voxel by the a priori 
                 mean time. To compute this parameter, we need to evaluate the squares of all
                 pair-wise cosines between the gradient directions*/
                m_SquaredCosines.SetSize(N,N);
                for( unsigned int r=0; r<N; ++r ){
                    for( unsigned int c=0; c<N; ++c ){
                        m_SquaredCosines[r][c]  = m_NeighboringDirections[r] * m_NeighboringDirections[c];
                        m_SquaredCosines[r][c] *= m_SquaredCosines[r][c];
                    }
                }
                /** We also have to precompute the sharpening
                 factors for the ODF and store the in an LUT */
                finsler_auxiliars::GenerateDeconvolutionFactors( m_LUT, m_bD, m_L );
                m_FRT = m_B; // In this case the input is directly the cOPDT 
                break;
            case TEST:
            case ISOTROPIC:
            case OPDFINV:
            case OPDFLOG:
            case DIFFUSIVITY:
                shmaths::computeSHEigMatrixSymmetric( m_L, eigenMatrix );
                m_FRT = m_B*m_FRT*eigenMatrix;
                break;
            default:
                itkExceptionMacro( << "Unsupported local cost computation method " << m_LocalCostType );
                break;
        }
        // The matrix m_B is the one required to compute the costs for each direction from
        // the SH coefficients at each voxel. Delete the memory previously allocated:
        delete[] theta;
        delete[] phi;
    }
    
    template< class TInputImage, class TOutputImage >
    void ComputeLocalCostFromSH< TInputImage, TOutputImage >
    ::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType )
    {
        // Iterators:
        ImageRegionConstIterator<InputImageType>        bit;  // Iterator for the input image
        ImageRegionIterator<OutputImageType>            it;   // Iterator for the output image
        // Input and output
        InputImageConstPointer   input   =  this->GetInput();
        OutputImagePointer       output  =  this->GetOutput();
        
        bit = ImageRegionConstIterator<InputImageType>( input,  outputRegionForThread );
        it  = ImageRegionIterator<OutputImageType>(     output, outputRegionForThread );
        
        // Precompute a couple of parameters of interest:
        unsigned int N = m_NeighboringDirections.size();
        unsigned int M = this->GetInput()->GetVectorLength();
        // Normalization in case we use logarithmic costs:
        const double log_norm = ( 12.566370614359172f / static_cast<double>(N) );
        // Auxiliar value to perform matrix products:
        MatrixType pixel(M,1);
        OutputPixelType op(N);
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
#ifdef CODE_FOR_DEBUG_EXPERIMENTAL
        // DEBUG
        typedef itk::Image<double,3>                     DebugImageType;
        typedef itk::ImageFileWriter<DebugImageType>     DebugWriterType;
        typedef itk::ImageRegionIterator<DebugImageType> DebugIteratorType;
        //...
        DebugImageType::Pointer debugIm = DebugImageType::New();
        debugIm->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
        debugIm->Allocate();
        debugIm->FillBuffer(0);
        debugIm->SetOrigin(this->GetInput()->GetOrigin());
        debugIm->SetSpacing(this->GetInput()->GetSpacing());
        debugIm->SetDirection(this->GetInput()->GetDirection());
        DebugIteratorType debugIt( debugIm, outputRegionForThread );
        debugIt.GoToBegin();
        //...
        DebugImageType::Pointer debugIm2 = DebugImageType::New();
        debugIm2->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
        debugIm2->Allocate();
        debugIm2->FillBuffer(0);
        debugIm2->SetOrigin(this->GetInput()->GetOrigin());
        debugIm2->SetSpacing(this->GetInput()->GetSpacing());
        debugIm2->SetDirection(this->GetInput()->GetDirection());
        DebugIteratorType debugIt2( debugIm2, outputRegionForThread );
        debugIt2.GoToBegin();
#endif
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
        for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it ){
            // Get the pixel and convert to a suitable format:
            InputPixelType  ip = bit.Get();
            for( unsigned int m=0; m<M; ++m )
                pixel[m][0] = ip[m];
            
            MatrixType att;
            MatrixType nor;
            double mean_diffusivity;
            double tau;
            
            switch( m_LocalCostType ){
                case ATTODF:
                   /* 
                    // OLD IMPLEMENTATION, CLONING JOHN MELONAKO'S CODE
                    // Convert to costs:
                    att = m_B*pixel;
                    nor = m_FRT*pixel;
                    // Compute the final output and set:
                    for( unsigned int n=0; n<N; ++n ){
                        double num   = ( att[n][0]>0 ? att[n][0] :        0       );
                        double den   = ( nor[n][0]>0 ? nor[n][0] : 0.01*att[n][0] );
                        double value = ( num + 0.00001f ) / ( den + 0.00001f );
                        // The 1000.0f x scaling is to avoid visualization issues:
                        op[n]        = 1000.0f * (value*value*value);
			//			op[n]        = 1000.0f * value;
                    }
                    break;
		   */
		    	
                    mean_diffusivity = pixel[0][0];
                    // Convert to costs:
                    nor = m_FRT*pixel;
                    // Compute the final output and set:
                    for( unsigned int n=0; n<N; ++n ){
                        double num   = 1.0f;
                        double den   = ( nor[n][0]>0 ? nor[n][0] : 0.01 );
                        double value = ( num + 0.00001f ) / ( den + 0.00001f );
                        // This is the difference with the previous case: we divide
                        // the cost by the mean diffusivity, so that voxels with high
                        // diffusivity propagate faster:
                        op[n]        = value*mean_diffusivity;//(value*value*value)/mean_diffusivity;
                    }
                    break;
             	    

                case OPDFINV:
                    // NEW IMPLEMENTATION BASED ON THE OPDF
                    // Convert to costs:
                    nor = m_FRT*pixel;
                    // Compute the final output and set:
                    for( unsigned int n=0; n<N; ++n ){
                        nor[n][0]   *= -0.025330295910584f;  // -1/4/pi^2
                        nor[n][0]   += 0.079577471545948;    // 1/4/pi
                        double num   = 1.0f;
                        double den   = ( nor[n][0]>0 ? nor[n][0] : 0.01 );
                        double value = ( num + 0.00001f ) / ( den + 0.00001f );
                        // Final value:
                        op[n]        = value;
						//op[n]        = (value*value*value);
                    }
                    break;
                case DIFFUSIVITY:
                    // OTHER NEW IMPLEMENTATION: we use the sharp shape of the
                    // ODF (see previous case) but multiply by the mean
                    // diffusivity to preserve the statistical meaning:
                    mean_diffusivity = pixel[0][0];
                    // Convert to costs:
                    nor = m_FRT*pixel;
                    // Compute the final output and set:
                    for( unsigned int n=0; n<N; ++n ){
                        nor[n][0]   *= -0.025330295910584f;  // -1/4/pi^2
                        nor[n][0]   += 0.079577471545948;    // 1/4/pi
                        double num   = 1.0f;
                        double den   = ( nor[n][0]>0 ? nor[n][0] : 0.01 );
                        double value = ( num + 0.00001f ) / ( den + 0.00001f );
                        // This is the difference with the previous case: we divide
                        // the cost by the mean diffusivity, so that voxels with high
                        // diffusivity propagate faster:
                        op[n]        = value*mean_diffusivity;//(value*value*value)/mean_diffusivity;
                    }
                    break;
                case MEANTIME:
                    // Preserve the basic diffusivity computed by
                    // itk::ComputeSHCoefficientsFilter:
                    mean_diffusivity = pixel[0][0];
                    // Fix the mean value of the ODF so that it has integral 1:
                    pixel[0][0] = 0.079577471545948; //1/(4*pi)
                    //---------------------------------------------------------------------------
                    // SHARPEN THE RESULTING ODF TAKING INTO ACCOUNT THE DIFFUSIVITY VALUE:
                    tau = this->SharpenODF( pixel, mean_diffusivity );
                    //---------------------------------------------------------------------------
                    // Compute the value of the cOPDT for each direction:
                    nor = m_FRT*pixel;
                    // Compute the final output (cost) and set:
                    for( unsigned int n=0; n<N; ++n ){
                        double num   = 1.0f;
                        double den   = ( nor[n][0]>0 ? nor[n][0] : 0.001 );
                        double value = ( num + 0.000001f ) / ( den + 0.000001f );
                        // Final value:
                        op[n]        = value;
                    }
                    break;
                    /*
                    // Sharpen the OPDF taking into account the diffusivity value. NOTE:
                    // since both the FRT and the L matrix are diagonal, they all commute
                    // with the sharpening matrix. This is the reason why we can sharp the
                    // coefficients of the attenuation signal, then compute the OPDF, and
                    // have the same result as firslty computing the OPDF, then sharpening
                    if( mean_diffusivity>0 ){
                        // Those regions not fulfilling our model (for example, those
                        // corresponding to the CSF) are encoded with a negative diffusivity
                        // value in itk::ComputeSHCoefficientsFilter
                        this->SharpenODF( pixel, mean_diffusivity );
                    }
                    else{
                        // In this case, we have detected a region very likely to correspond
                        // to the CSF, or at least a region that no longer fulfills our
                        // model. Hence, we consider here an isotropic behavior:
                        pixel.fill(0);
                        // And correct the value of the diffusivity
                        mean_diffusivity = 1;
                    }
                    // Convert to OPDF:
                    nor = m_FRT*pixel;
                    // Compute the final output and set:
                    for( unsigned int n=0; n<N; ++n ){
                        nor[n][0]   *= -0.025330295910584f;
                        nor[n][0]   += 0.079577471545948;
                        double num   = 1.0f;
                        double den   = ( nor[n][0]>0 ? nor[n][0] : 0.01 );
                        nor[n][0]    = den;
                        double value = ( num + 0.00001f ) / ( den + 0.00001f );
                        // Inverse of the sharpened OPDF:
                        op[n]        = value;
                    }
                    // Compute the a priori mean diffusion time; note we don't need
                    // the OPDF to be normalized to sum to 1, since any constant 
                    // factor multiplying it won't effect the value of tau (we have
                    // the OPDF both in the numerator and the denominator of tau).
                    tau = this->ComputeAPrioriMeanTime( nor );
                    //-------------------------------------------------------------------------------------------------
                    //-------------------------------------------------------------------------------------------------
                    //-------------------------------------------------------------------------------------------------
#ifdef CODE_FOR_DEBUG_EXPERIMENTAL
                    // DEBUG
                    debugIt.Set( tau );
                    ++debugIt;
                    debugIt2.Set( mean_diffusivity );
                    ++debugIt2;
                    mean_diffusivity = itk::NumericTraits<double>::One;
#endif
                    //-------------------------------------------------------------------------------------------------
                    //-------------------------------------------------------------------------------------------------
                    //-------------------------------------------------------------------------------------------------
                    // Assign the final value to the pixel:
                    tau = vcl_sqrt( tau/mean_diffusivity );
                    for( unsigned int n=0; n<N; ++n )
                        op[n] *= tau;
                    break;
                     */
                    
                case OPDFLOG:
                    // Convert to costs:
                    nor = m_FRT*pixel;
                    // Compute the final output and set:
                    for( unsigned int n=0; n<N; ++n ){
                        nor[n][0]   *= -0.025330295910584f;
                        nor[n][0]   += 0.079577471545948;
                        double val   = ( nor[n][0]>0.0001 ? nor[n][0] : 0.0001 );
                        val         *= log_norm;
                        val          = ( val<0.9999 ? val : 0.9999 );
                        op[n]        = -3.0f * vcl_log( val );
                    }
                    break;
                case ISOTROPIC:
                    for( unsigned int n=0; n<N; ++n ){
                        op[n] = 1.0f;
                    }
                    break;
                case TEST:
                    // DUMB TEST FOR DEBUG
                    if( bit.GetIndex()[0]>2 ){
                        for( unsigned int n=0; n<N; ++n ){
                            DirectionType d1 = m_NeighboringDirections[n];
                            d1.Normalize();
                            DirectionType d2;
                            d2[0] = 1.0f;
                            d2[1] = 0.0f;
                            d2[2] = 0.0f;
                            double temp = d1*d2;
                            temp = ( temp>itk::NumericTraits<double>::Zero ? temp : -temp );
                            if( temp>0.9 && temp<1.1 )
                                op[n] = 1.0f;
                            else
                                op[n] = 200.0f;
                        }
                    }
                    else{
                        for( unsigned int n=0; n<N; ++n )
                            op[n] = 5.0f;
                    }
                    break;
                default:
                    itkExceptionMacro( << "Unsupported type of local cost");
                    break;
            }
            it.Set( op );
        }
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
#ifdef CODE_FOR_DEBUG_EXPERIMENTAL
        // DEBUG
        DebugWriterType::Pointer debugWriter = DebugWriterType::New();
        debugWriter->SetInput( debugIm );
        debugWriter->SetFileName("/Users/atriveg/Downloads/debugPriorTau.nrrd");
        debugWriter->Update();
        //...
        DebugWriterType::Pointer debugWriter2 = DebugWriterType::New();
        debugWriter2->SetInput( debugIm2 );
        debugWriter2->SetFileName("/Users/atriveg/Downloads/debugD.nrrd");
        debugWriter2->Update();
        std::cout << "Finished prior computations" << std::endl;
#endif
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------
    }
    
    template< class TInputImage, class TOutputImage >
    double ComputeLocalCostFromSH< TInputImage, TOutputImage >
    ::ComputeAPrioriMeanTime( const MatrixType& opdf )
    {
        /** This function is used only in the "MEANTIME" case. From a vector representing
         the ODF, it computes the mean of the (exponential) prior distribution of the
         diffusion time tau according to the model we use. Note this mean time is normalized
         by the a priori mean time for free diffusion, so the value returned cannot be
         considered a physical parameter (it has no units). */
        
        // The number of diffusion directions considered:
        unsigned int N = m_NeighboringDirections.size();
        // The value to be returned:
        double tau = itk::NumericTraits<double>::Zero;
        
        for( unsigned int i=0; i<N; ++i ){
            // For each direction considered, compute the mean contribution 
            // of the remaining directions to the current one:
            double  mean_cosines = itk::NumericTraits<double>::Zero;
            for( unsigned int j=0; j<N; ++j )
                mean_cosines += ( opdf[j][0] * m_SquaredCosines[i][j] );
            // The mean time tau is the average of the inverse of these contributions:
            tau += (   ( opdf[i][0] + 1.0e-15f )  /  ( mean_cosines + 1.0e-15f )   );
        }
        
        return(tau);
    }
    
    template< class TInputImage, class TOutputImage >
    double ComputeLocalCostFromSH< TInputImage, TOutputImage >
    ::SharpenODF( MatrixType& shcoeffs, const double& attenuation )
    {
        /** This function is based on the following paper:
         A. Tristan-Vega, S. Aja-Fernandez, and C.-F. Westin
         "Deblurring of Probabilistic ODFs in Quantitative Diffusion MRI"
         In Proceedings of the IEEE Intl. Sym. on Biomed. Im., ISBI 2012
         Barcelona (Spain), pp. 932-935
         */
        // We start by computing the basic diffusivity bD from the mean
        // attenuation. This is done by inverting eq. (6) in the
        // aforementioned paper:
        //    attenuation = sqrt(pi/4) * erf( sqrt(bD) ) / sqrt(bD)
        // To do so, we have fitted f(t) = sqrt(pi/4) * erf(t) / t
        // for different intervals using series expansions, so that
        // the relative error |f(t)-fit|/|f(t)| is AT MOST <0.5%:
        //
        //    t in (0,1/2) => f(t) in (0.922562,1):
        //         f(t) ~= 1 - t^2/3 + t^4/10 (taylor series)
        //    t in (1/2,1) => f(t) in (0.746824,0.922562):
        //         f(t) ~= -0.0900(t-3/4)^2 - 0.3551(t-3/4) + 0.8403 (empirical)
        //    t in (1,3/2) => f(t) in (0.570792,0.746824):
        //         f(t) ~=  0.0723(t-5/4)^2 - 0.3535(t-5/4) + 0.6543 (empirical)
        //    t in (3/2,2) => f(t) in (0.441040,0.570792):
        //         f(t) ~=  0.1002(t-7/4)^2 - 0.2592(t-7/4) + 0.4997 (empirical)
        //    t>2          => f(t) < 0.441040:
        //         f(t) ~= ~= sqrt(pi/4) / t  (series expansion)
        
        // We have to invert each branch for bD
        double bD = itk::NumericTraits<double>::Zero;
        if( attenuation > 0.922562f ){
            // t^2 = 5/3 - sqrt( (5/3)^2 - 10(1-att) )
            bD  = -10*(1-attenuation);
            bD += 2.777777777777778;
            bD  = -vcl_sqrt( bD );
            bD += 1.666666666666667;
            bD  = vcl_sqrt(bD);
        }
        else if( attenuation > 0.746824f ){
            // att = a*s^2 + b*s + c => a*s^2 + b*s + (c-att) = 0; (s=t-t0)
            //   t = t0 + ( -b +- sqrt( b^2 - 4a(c-att) ) ) / 2a
            const double t0 =  0.75f;
            const double a  = -0.089962815105939f;
            const double b  = -0.355142505587604f;
            const double c  =  0.840334139440334f;
            // Must choose the minus sign so that the result is always positive:
            bD = t0 + ( -b - vcl_sqrt(b*b-4*a*(c-attenuation)) )/a/2;
        }
        else if( attenuation > 0.570792f ){
            // att = a*s^2 + b*s + c => a*s^2 + b*s + (c-att) = 0; (s=t-t0)
            //   t = t0 + ( -b +- sqrt( b^2 - 4a(c-att) ) ) / 2a
            const double t0 =  1.25f;
            const double a  =  0.072268143989957f;
            const double b  = -0.353545105875437f;
            const double c  =  0.654336220350557f;
            // Must choose the minus sign so that the result is always positive:
            bD = t0 + ( -b - vcl_sqrt(b*b-4*a*(c-attenuation)) )/a/2;
        }
        else if( attenuation > 0.441040f ){
            // att = a*s^2 + b*s + c => a*s^2 + b*s + (c-att) = 0; (s=t-t0)
            //   t = t0 + ( -b +- sqrt( b^2 - 4a(c-att) ) ) / 2a
            const double t0 =  1.75f;
            const double a  =  0.100169288804115f;
            const double b  = -0.259230410348181f;
            const double c  =  0.499671645686046f;
            // Must choose the minus sign so that the result is always positive:
            bD = t0 + ( -b - vcl_sqrt(b*b-4*a*(c-attenuation)) )/a/2;
        }
        else
            bD = 0.886226925452758f/attenuation; // sqrt(pi/4) = 0.886226925452758
        // We have computed t = sqrt(bD), so we have to compute the
        // square of this magnitude:
        bD *= bD;
        if( bD<5 )
            bD = 5;
        /** ********************************************************************** */
        // AT THIS POINT, WE HAVE A VALUE OF bD SUCH THAT:
        //   attenuation = sqrt(pi/4) * erf( sqrt(bD) ) / sqrt(bD);
        // with a relative error |bD_estimated-bD|/|bD| at most <1% (this has
        // been empirically checked with a matlan script test_fit.m)
        /** ********************************************************************** */
        // Now we have computed the actual diffusivity, and will be able to find
        // the linear correction coefficients for the SH; How many entries do
        // we have in the LUT?
        const unsigned int res = m_bD.size();
        // We use a nearest-neighbor approach to find the position in the LUT:
        unsigned int pos = 0;
        if( bD<=m_bD[0] )          // The diffusivity is too small!
            pos = 0;
        else if( bD>=m_bD[res-1] ) // The diffusivity is too large!
            pos = res-1;
        else{
            for( pos=0; pos<res; ++pos ){
                if( m_bD[pos]>bD )
                    break;
            }
        }
        // Get the appropriate correction factors:
        std::vector<double> factors = m_LUT[pos];
        // Apply to the SH coefficients to get the sharp OPDF:
        // BEWARE: For the sake of performance, we do not check the size of
        //         shcoeffs, which we assume has the proper size for the
        //         following loop:
        pos = 1; // Skip the DC component, which is always sqrt(1/4/PI)
        for( unsigned int l=2; l<=m_L; l+=2 ){ // Skip the DC component, which is always sqrt(1/4/PI)
            for( int m=-(int)l; m<=(int)l; ++m )
                shcoeffs( pos++, 0 ) *= factors[l/2];
        }
        // Done
        return bD;
    }
    
} // end namespace itk


#endif
