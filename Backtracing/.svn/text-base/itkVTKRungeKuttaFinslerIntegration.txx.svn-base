/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkVTKRungeKuttaFinslerIntegration.txx,v $
 Language:  C++
 Date:      $Date: 2006/01/11 19:43:31 $
 Version:   $Revision: 1.21 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef _itkVTKRungeKuttaFinslerIntegration_txx
#define _itkVTKRungeKuttaFinslerIntegration_txx

#include "itkVTKRungeKuttaFinslerIntegration.h"

namespace itk
{
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::VTKRungeKuttaFinslerIntegration()
    {
        m_MinimumLength = 10.0f;
        m_MaximumLength = 800.0f;
        m_StoppingCurvature = 0.8f;
        m_StepLength = 0.5f;
        m_Label = 1;
        
        m_Interpolate       = NULL;
        m_VectorInterpolate = NULL;
        m_Streamlines       = NULL;
        
        m_NormalizedCosts   = NULL;
        
        this->SetNumberOfRequiredInputs(3);
        // A nice further improvement would be to make the filter multi-threaded:
        this->SetNumberOfThreads(1);       
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    void VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::SetArrivals(const ArrivalsType* arrivals)
    {
        this->Superclass::SetNthInput( 0, const_cast<ArrivalsType*>( arrivals ) );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    void VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::SetCosts(const CostType* costs)
    {
        this->Superclass::SetNthInput( 1, const_cast<CostType*>( costs ) );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    void VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::SetLabels(const LabelType* labels)
    {
        this->Superclass::SetNthInput( 2, const_cast<LabelType*>( labels ) );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    const typename VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::ArrivalsType*
    VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::GetArrivals()
    {
        return static_cast<ArrivalsType*> ( this->Superclass::GetInput(0) );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    const typename VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::CostType*
    VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::GetCosts()
    {
        return static_cast<CostType*> ( this->Superclass::GetInput(1) );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    const typename VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::LabelType*
    VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::GetLabels()
    {
        return static_cast<LabelType*> ( this->Superclass::GetInput(2) );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    void VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::Update()
    {
        // PREPARATION AND SETTING-UP ----------------------------------------------------
        this->Superclass::Update();
        
        const ArrivalsType* arrivals = this->GetArrivals();
        const CostType*     costs    = this->GetCosts();
        const LabelType*    labels   = this->GetLabels();
        
        if( !arrivals )
            itkExceptionMacro( << "Please, set the optimal arrival directions map" );
        if( !costs )
            itkExceptionMacro( << "Please, set the costs map" );
        if( !labels )
            itkExceptionMacro( << "Please, set labels map" );
        if( arrivals->GetSource() )
            arrivals->GetSource()->Update();
        if( costs->GetSource() )
            costs->GetSource()->Update();
        if( labels->GetSource() )
            labels->GetSource()->Update();
        
        // All the set of streamlines: I create a new set with each new call to
        // Update(); since I use smart pointers, the memory associated to the old
        // pointers should be erased as needed:
        m_Streamlines = vtkSmartPointer<vtkPolyData>::New();
        // All the points in all streamlines:
        vtkSmartPointer<vtkPoints>     points      = vtkSmartPointer<vtkPoints>::New();
        // The set of scalars to color the fibers:
        vtkSmartPointer<vtkFloatArray> scalars     = vtkSmartPointer<vtkFloatArray>::New();
        // The set of indices of points in the current streamline:
        vtkSmartPointer<vtkIdList>     current     = vtkSmartPointer<vtkIdList>::New();
        
        m_Streamlines->Initialize();
        m_Streamlines->Allocate();
        
        points->Initialize();
        points->Allocate(0);
        
        scalars->Initialize();
        scalars->Allocate(0);
        
        if( !m_Interpolate )
            m_Interpolate       = InterpolateType::New();
        if( !m_VectorInterpolate )
            m_VectorInterpolate = VectorInterpolateType::New(); 
        m_Interpolate->SetInputImage( costs );
        m_VectorInterpolate->SetInputImage( arrivals );
        // END PREPARATION AND SETTING-UP ------------------------------------------------
        
        // ACTUALLY DO THE JOB!
        unsigned long numberOfPoints = 0; // Global number of points in the polydata
        IteratorType it( labels, labels->GetLargestPossibleRegion() );
        for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
            // For each point in the labelmap image, we check if this is
            // a target to trace from.
            if( it.Get() != m_Label )
                continue;
            // If we reach this point, this is a target point. We need
            // its physical position:
            PointType currentTarget;
            for( unsigned int k=0; k<TCostMapType::ImageDimension; ++k )
                currentTarget[k] = it.GetIndex()[k] * labels->GetSpacing()[k];
            currentTarget  = labels->GetDirection() * currentTarget;
            for( unsigned int k=0; k<TCostMapType::ImageDimension; ++k )
                currentTarget[k] += labels->GetOrigin()[k];
            // Now currentTarget is a physical point and we can start the tracking:
            unsigned int causeOfStop = 0;
            double length = this->StreamNewFiberBundle( currentTarget, points, scalars, 
                                                            current, numberOfPoints, causeOfStop );
            // The new streamline will be included only if it fulfills all the requirements:
            if( length>=m_MinimumLength && length<=m_MaximumLength ){
                m_Streamlines->InsertNextCell( VTK_POLY_LINE, current );
                // In  case we must write the normalized costs, do it now. It's not
                // really elegant using "SetPixel()" instead of using an iterator.
                // However, since the ROI is assumed to be much smaller than the image,
                // and since we have to compute a whole fiber for each step of the
                // iterator, the computational overload is completely negligible.
                if( m_NormalizedCosts ){
                    // The cost is only trusted if the stopping criterion is considered
                    // "normal", i.e., we are at the seeding region. This implies the
                    // causeOfStop is either 4 or 5; otherwise, a -1 normalized cost
                    // is placed:
                    typename CostType::IndexType currentIndex;
                    if( m_NormalizedCosts->TransformPhysicalPointToIndex( currentTarget, currentIndex ) ){
                        if( (causeOfStop==4) || (causeOfStop==5) )
                            m_NormalizedCosts->SetPixel( currentIndex, 10000 * it.Get()/length );
                        else
                            m_NormalizedCosts->SetPixel( currentIndex, -itk::NumericTraits<typename CostType::PixelType>::One );
                    }
                }
            }
            else{
                // The streamline didn't pass the test. A -1 cost should be
                // placed in the output costs map if needed:
                if( m_NormalizedCosts ){
                    typename CostType::IndexType currentIndex;
                    if( m_NormalizedCosts->TransformPhysicalPointToIndex( currentTarget, currentIndex ) ){
                        m_NormalizedCosts->SetPixel( currentIndex, -itk::NumericTraits<typename CostType::PixelType>::One );
                    }
                }
            }
        }
        // Once all target points have been tracked, we can tell the polydata
        // what points it has to use:
        m_Streamlines->SetPoints( points );
        m_Streamlines->GetPointData()->SetScalars( scalars );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    double VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::StreamNewFiberBundle(
        typename VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType,TCostMapType,TLabelImageType>::PointType currentTarget,
        vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkFloatArray> scalars, 
        vtkSmartPointer<vtkIdList> current, unsigned long& numberOfPoints, unsigned int& causeOfStop )
    {
        // Parameters of interest:
        const double finsler_cost_tol    = ( m_Interpolate->IsInsideBuffer(currentTarget)
                                                    ?  1.0e-3 * m_Interpolate->Evaluate(currentTarget)
                                                    :  1.0e-3       );
        const double minimum_jump_length = 1.0e-2*m_StepLength;
        double length    = 0.0f;
        double curvature = 2*m_StoppingCurvature;
        bool   stop      = false;
        // Points for the integration (We initialize them to the current
        // target to avoid compiler issues, but their value will not be
        // used in the first iterations):
        PointType nextTarget     = currentTarget;
        PointType prevTarget     = currentTarget;
        // This variable contains the numbe of points in the current streamline:
        unsigned long pointsInThisStream = 0;
        // The current streamline must be initialized:
        current->Initialize();
        current->Allocate(0);
        // Now we can start the loop to actually compute the streamline:
        do{
            causeOfStop = 0; // Auxiliar to decide wether to stop
            // Include the current target point in the current streamline; though it
            // should never happen, we perform a test to check that this target is
            // inside the buffered regions (otherwise, this is a bad taget and the
            // streamline will be discarded
            if( m_Interpolate->IsInsideBuffer(currentTarget) && m_VectorInterpolate->IsInsideBuffer(currentTarget) ){
                // The point is OK, so we can insert it in the polydata, its cost in
                // the scalar list (to define its color), and its Id in the list for
                // the current streamline.
                //---
                // We will write the fiber bundles to the disk always in RAS coordinates. We do
                // so because the vtk formats used to store the fiber bundles do not encode the
                // space (RAS/LPS...), and Slicer will assume the fibers are always in the RAS
                // space. However, ITK works always in the left-posterior-superior (LPS) space,
                // so we need to convert from LPS to RAS. Fortunately, this implies only shifting
                // the 'x' and 'y0 components of the 3D coordinate:
                points->InsertNextPoint( -currentTarget[0], -currentTarget[1], currentTarget[2] );
                scalars->InsertNextValue( (float)(m_Interpolate->Evaluate(currentTarget) ) );
                current->InsertNextId( numberOfPoints++ ); // One more point!
                // And update the length of the fiber:
                length += ( currentTarget - prevTarget ).GetNorm();
                pointsInThisStream++;
                // Now, use the Runge-Kutta method implementation to integrate the next point:
                curvature = ComputeNextPointInStream( prevTarget, currentTarget,
                                                        nextTarget, pointsInThisStream );
            }
            else
                causeOfStop = 1;
            // Now, check if this point is correct; otherwise, we will return. We have
            // to exit if, for some reason, the point is outside of the buffered region
            // of any of the images:
            if( causeOfStop==0 )
                causeOfStop = ( m_Interpolate->IsInsideBuffer( nextTarget ) ? 0 : 2 );
            if( causeOfStop==0 )
                causeOfStop = ( m_VectorInterpolate->IsInsideBuffer( nextTarget ) ? 0 : 3 );
            // We also have to stop if the cost is such that we are in the seeding
            // region or outside of the mask (if used):
            if( causeOfStop==0 )
                causeOfStop = ( m_Interpolate->Evaluate(nextTarget) >= finsler_cost_tol ? 0 : 4 );
            // Stop if we have repeated points:
            if( causeOfStop==0 )
                causeOfStop = ( (nextTarget-currentTarget).GetNorm()>minimum_jump_length ? 0 : 5 );
            // Stop if the lenhgth of the fiber is excessive:
            if( causeOfStop==0 )
                causeOfStop = ( length<m_MaximumLength ? 0 : 6 );
            // Finally, we have to exit if the curvature radius is too small:
            if( causeOfStop==0 )
                causeOfStop = ( curvature>=m_StoppingCurvature ? 0 : 7 );
            // This is an "emergency stop" to avoid infinite loops:
            if( causeOfStop==0 )
                causeOfStop = ( pointsInThisStream <= MAXIMUM_NUMBER_OF_POINTS_PER_STREAM ? 0 : 8 );
            // Gather all possible causes of stopping:
            if( causeOfStop != 0 )
                stop = true;
            else{
                // Slide the window to keep on integrating:
                prevTarget     = currentTarget;
                currentTarget  = nextTarget;
            }
        }
        while( !stop );
          
        return( length );
    }
    
    template <class TArrivalDirectionsType, class TCostMapType, class TLabelImageType>
    double VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType, TCostMapType, TLabelImageType>
    ::ComputeNextPointInStream(
        const typename VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType,TCostMapType,TLabelImageType>::PointType& prevTarget,
        const typename VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType,TCostMapType,TLabelImageType>::PointType& currentTarget,
        typename VTKRungeKuttaFinslerIntegration<TArrivalDirectionsType,TCostMapType,TLabelImageType>::PointType& nextTarget,
        const unsigned long pointsInThisStream )
    {
        // First, compute the initial trajectory for Euler integration. For the first
        // point we take the vector in the directions map corresponding to the heading
        // point in the streamline. Otherwise, this method will compute the value and
        // subsequently store it in "trajectories".
        StepType step;
        if( m_VectorInterpolate->IsInsideBuffer( currentTarget ) ){
            typename VectorInterpolateType::OutputType temp = m_VectorInterpolate->Evaluate( currentTarget );
            for( unsigned int d=0; d<TCostMapType::ImageDimension; ++d )
                step[d] = temp[d];
        }
        else{
            // This situation should never arise. We return with a
            // wrong curvature to tell the algorithm to stop:
            return( m_StoppingCurvature/2 );
        }
        if( step.GetNorm()<1e-3 ){
            // This situation should never arise. We return with a
            // wrong curvature to tell the algorithm to stop:
            return( m_StoppingCurvature/2 );
        }
        else
            step.Normalize();
        // Note the trajectory is assumed to be already in real world units,
        // so we only have to move m_StepLength units in the direction
        // provided. This is the Euler integration step:
        nextTarget = currentTarget + ( step * m_StepLength );
        // Now we are gonna retrieve the displacement vector in the udpated
        // position; in case that position is outside the buffered region,
        // we keep the corrent step:
        StepType step_next;
        if( m_VectorInterpolate->IsInsideBuffer( nextTarget ) ){
            typename VectorInterpolateType::OutputType temp = m_VectorInterpolate->Evaluate( nextTarget );
            for( unsigned int d=0; d<TCostMapType::ImageDimension; ++d )
                step_next[d] = temp[d];
        }
        else
            step_next = step;
        if( step_next.GetNorm()<1e-3 ){
            // We have reached the seeding region
            step_next = step;
        }
        else
            step_next.Normalize();
        // Find the definitive position of the next point with the
        // resulting step:
        nextTarget = currentTarget + ( (step+step_next) * (0.5f*m_StepLength) );
        // Finally, compute the current curvature that we will return
        // to be used as a stop criterion:
        double curvature = 2.0f * m_StoppingCurvature;
        if( pointsInThisStream > 1 ){
            // We have enough points to compute the curvature and we do so.
            // Otherwise, we keep on integratig until we have at least 2 points.
            StepType kv2 = prevTarget - currentTarget;
            StepType kv1 = currentTarget - nextTarget;
            
            double   kl2 = kv2.GetNorm();
            double   kl1 = kv1.GetNorm();
            
            kv2.Normalize();
            kv1.Normalize();
            
            // The curvature is:
            StepType kn  = ( kv2 - kv1 ) * ( 2.0f/(kl1+kl2) );
            curvature = kn.GetNorm();
            
            // And the radius is its inverse:
            if( curvature > 1e-3 )
                curvature = 1.0f/curvature;
            else
                curvature = 2.0f * m_StoppingCurvature;
        }
        return( curvature );
    }
    
} // end namespace itk


#endif
