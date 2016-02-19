/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkVTKRungeKuttaFinslerIntegration.h,v $
 Language:  C++
 Date:      $Date: 2006/03/27 17:01:10 $
 Version:   $Revision: 1.15 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef __itkVTKRungeKuttaFinslerIntegration_h
#define __itkVTKRungeKuttaFinslerIntegration_h

#include "itkProcessObject.h"
#include "itkImage.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVector.h"
#include "itkFixedArray.h"
#include <vector>

#define MAXIMUM_NUMBER_OF_POINTS_PER_STREAM 100000

namespace itk
{
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    class ITK_EXPORT VTKRungeKuttaFinslerIntegration : public ProcessObject
    {
    public:
        /** Standard class typedefs. */
        typedef VTKRungeKuttaFinslerIntegration           Self;
        typedef SmartPointer<Self>                        Pointer;
        typedef SmartPointer<const Self>                  ConstPointer;
        typedef ProcessObject                             Superclass;
                
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro( VTKRungeKuttaFinslerIntegration, ProcessObject );
        
        /** Convenient typedefs for simplifying declarations. */
        typedef TArrivalDirectionsType                    ArrivalsType;
        typedef typename ArrivalsType::Pointer            ArrivalsPointer;
        typedef typename ArrivalsType::ConstPointer       ArrivalsConstPointer;
        
        typedef TCostMapType                              CostType;
        typedef typename CostType::Pointer                CostPointer;
        typedef typename CostType::ConstPointer           CostConstPointer;
        
        typedef TLabelImageType                           LabelType;
        typedef typename LabelType::Pointer               LabelPointer;
        typedef typename LabelType::ConstPointer          LabelConstPointer;
        
        typedef typename ArrivalsType::PixelType          VectorType;
        typedef typename CostType::PixelType              TimeType;
        typedef typename LabelType::PixelType             TargetType;
        
        typedef typename ArrivalsType::RegionType         RegionType;
        typedef typename ArrivalsType::SizeType           SizeType;
        typedef typename ArrivalsType::IndexType          IndexType;
        typedef typename ArrivalsType::PointType          PointType;
        typedef typename ArrivalsType::SpacingType        SpacingType;
        typedef typename ArrivalsType::DirectionType      DirectionType;
        
        typedef itk::Vector<double,TCostMapType::ImageDimension>    StepType;
        
        /** These interpolators will be needed to track the finbers. */
        typedef itk::LinearInterpolateImageFunction<CostType>       InterpolateType;
        typedef typename InterpolateType::Pointer                   InterpolatePointer;
        
        typedef itk::VectorLinearInterpolateImageFunction<ArrivalsType>   VectorInterpolateType;
        typedef typename VectorInterpolateType::Pointer                   VectorInterpolatePointer;
        
        /** Type to iterate through the map of potential targets. */
        typedef itk::ImageRegionConstIteratorWithIndex<LabelType>   IteratorType;
        
        vtkSmartPointer<vtkPolyData> GetStreamlines( void )
        {
            return m_Streamlines;
        }
        
        void SetArrivals( const ArrivalsType* );
        void SetCosts( const CostType* );
        void SetLabels( const LabelType* );
        
        const ArrivalsType* GetArrivals( void );
        const CostType* GetCosts( void );
        const LabelType* GetLabels( void );
        
        itkSetMacro( MinimumLength, double );
        itkGetMacro( MinimumLength, double );
        
        itkSetMacro( MaximumLength, double );
        itkGetMacro( MaximumLength, double );
        
        itkSetMacro( StoppingCurvature, double );
        itkGetMacro( StoppingCurvature, double );
        
        itkSetMacro( StepLength, double );
        itkGetMacro( StepLength, double );
        
        itkSetMacro( Label, unsigned int );
        itkGetMacro( Label, unsigned int );
        
        void SetNormalizedCosts( CostPointer& costs ){
            m_NormalizedCosts = costs;
            return;
        }
        
        CostPointer& GetNormalizedCosts( void ){
            return m_NormalizedCosts;
        }
        
        virtual void Update( void );
    protected:
        VTKRungeKuttaFinslerIntegration();
        virtual ~VTKRungeKuttaFinslerIntegration() {}
        // This method is the one that actually streams the fiber bundles
        double StreamNewFiberBundle( PointType, vtkSmartPointer<vtkPoints>, vtkSmartPointer<vtkFloatArray>,
                                                vtkSmartPointer<vtkIdList>, unsigned long&, unsigned int& );
        // This method computes one single step of the Runge-Kutta integration:
        double ComputeNextPointInStream( const PointType&, const PointType&, PointType&,
                                                 const unsigned long );
    private:
        VTKRungeKuttaFinslerIntegration(const Self&); // purposely not implemented
        void operator=(const Self&);                  // purposely not implemented
        // I'm using smart pointers to avoid worrying about memory management:
        vtkSmartPointer<vtkPolyData>                      m_Streamlines;
        
        // Interpolators:
        InterpolatePointer                                m_Interpolate;
        VectorInterpolatePointer                          m_VectorInterpolate;
        
        // Specific parameters:
        double m_MinimumLength;
        double m_MaximumLength;
        double m_StoppingCurvature;
        double m_StepLength;
        unsigned int m_Label;
        
        // This pointer is used to store the costs at the target ROI locations
        // normalized by the arc length of the corresponding fuber bundle. In
        // some situations it is possible that the backtracing algorithm does
        // not find a parth connecting the target to the seed due to an excessive
        // curvature or arc length. In these situations, the normalized cost is
        // set to "-1". Note the use of this output map is optional.
        // VERY IMPORTANT: No size/memory allocation checking is performed, so
        // the user is responsible of providing a pointer with the same size,
        // origin, spacing, and orientation as the input costs map used.
        CostPointer m_NormalizedCosts;
    };
    
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVTKRungeKuttaFinslerIntegration.txx"
#endif

#endif
