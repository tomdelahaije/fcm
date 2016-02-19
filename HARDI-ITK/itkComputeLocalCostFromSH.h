/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkComputeLocalCostFromSH.h,v $
 Language:  C++
 Date:      $Date: 2006/03/27 17:01:10 $
 Version:   $Revision: 1.15 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef __itkComputeLocalCostFromSH_h
#define __itkComputeLocalCostFromSH_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkArray2D.h"
#include "sphericalHarmonics.h"

namespace itk
{

    template <class TInputImage, class TOutputImage>
    class ITK_EXPORT ComputeLocalCostFromSH : public ImageToImageFilter< TInputImage, TOutputImage >
    {
    public:
        /** Standard class typedefs. */
        typedef ComputeLocalCostFromSH                Self;
        /** Convenient typedefs for simplifying declarations. */
        typedef TInputImage                           InputImageType;
        typedef typename InputImageType::Pointer      InputImagePointer;
        typedef typename InputImageType::ConstPointer InputImageConstPointer;
        typedef TOutputImage                          OutputImageType;
        typedef typename OutputImageType::Pointer     OutputImagePointer;
        
        /** Standard class typedefs. */
        typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
        typedef SmartPointer<Self>                                   Pointer;
        typedef SmartPointer<const Self>                             ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro( ComputeLocalCostFromSH, ImageToImageFilter );
        
        /** Image typedef support. */
        typedef typename InputImageType::PixelType           InputPixelType;
        typedef typename OutputImageType::PixelType          OutputPixelType;
        typedef typename InputImageType::RegionType          InputImageRegionType;
        typedef typename InputImageType::SizeType            InputImageSizeType;
        typedef typename InputImageType::IndexType           InputImageIndexType;
        typedef typename OutputImageType::RegionType         OutputImageRegionType;
        
        /** Types to store the directions of interest */
        typedef itk::CovariantVector<float,TInputImage::ImageDimension> DirectionType;
        
        void SetNeighboringDirections( std::vector<DirectionType> dirs )
        {
            m_NeighboringDirections = dirs;
        }
        
        /** How do we compute the local cost */
        typedef enum{ ATTODF, OPDFINV, DIFFUSIVITY, MEANTIME, OPDFLOG, ISOTROPIC, TEST } LocalCostType;
        
        void SetComputeATTOPDF(void)
        {
            m_LocalCostType = ATTODF;
        }
        
        void SetComputeOPDFINV(void)
        {
            m_LocalCostType = OPDFINV;
        }
        
        void SetComputeDIFFUSIVITY(void)
        {
            m_LocalCostType = DIFFUSIVITY;
        }
        
        void SetComputeMEANTIME(void)
        {
            m_LocalCostType = MEANTIME;
        }
        
        void SetComputeOPDFLOG(void)
        {
            m_LocalCostType = OPDFLOG;
        }

        void SetComputeISOTROPIC(void)
        {
            m_LocalCostType = ISOTROPIC;
        }
       
        void SetComputeTEST(void)
        {
            m_LocalCostType = TEST;
        }
        
        // Type for matrix computations (vnl matrix);
        typedef shmaths::SHMatrixType                        MatrixType;
        
        // Type to store the squares of the pair-wise cosines
        typedef itk::Array2D<float>                          SquaredCosinesType;
    protected:
        ComputeLocalCostFromSH();
        virtual ~ComputeLocalCostFromSH() {}
        void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType );
        void BeforeThreadedGenerateData( void );
        void GenerateOutputInformation( void );
        double ComputeAPrioriMeanTime( const MatrixType& );
        // This method is only used in the MEANTIME case (at this moment; it could be
        // used whenever the inverse of the cOPDT, or even the OPDT, is used). It
        // maniopulates the SH coefficients to deconvolve the OPDF and avoid the
        // blurring coming from the low-pass filtering inside the circle q=q_0:
        double SharpenODF( MatrixType&, const double& );
    private:
        ComputeLocalCostFromSH(const Self&);    // purposely not implemented
        void operator=(const Self&);            // purposely not implemented
        
        // The directions for which the cost will be computed:
        std::vector<DirectionType> m_NeighboringDirections;
        // The matrixes to compute the values from the SH coefficients:
        MatrixType                 m_B;   // SH matrix
        MatrixType                 m_FRT; // The Funk-Radon transform matrix
        // How do we compute the local cost
        LocalCostType              m_LocalCostType;
        // Store the squared cosines (only useful for the MEANTIME case)
        SquaredCosinesType         m_SquaredCosines;
        // These structures contain a look-up table that stores the
        // factors we have to apply to the SH coefficients of the OPDF
        // to avoid the blurring due to the low-pass filtering. In this
        // version of the software tehyr are used only in the MEANTIME
        // case, though they can be used whenever the cOPDT (or even the
        // OPDT) is involved.
        std::vector<double>                 m_bD;
        std::vector< std::vector<double> >  m_LUT;
        // The order of the SH:
        unsigned int m_L;
    };
    
    
    
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeLocalCostFromSH.txx"
#endif

#endif
