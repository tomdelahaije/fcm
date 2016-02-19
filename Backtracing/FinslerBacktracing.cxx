/*=========================================================================
 
 Program:   Diffusion Applications
 Language:  C++
 Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/dwiNoiseFilter/dwiNoiseFilter.cxx $
 Date:      $Date: 2008-11-25 18:46:58 +0100 (Tue, 25 Nov 2008) $
 Version:   $Revision: 7972 $
 
 Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.
 
 See License.txt or http://www.slicer.org/copyright/copyright.txt for details.
 
 ==========================================================================*/
#include <iostream>
#include <algorithm>
#include <string>
#include <itkMetaDataObject.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkVectorImage.h>

#ifdef _WIN32
// to pick up M_SQRT2 and other nice things...
#define _USE_MATH_DEFINES
#endif
#include <math.h>


#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include "itkCastImageFilter.h"
#include "FinslerBacktracingCLP.h"

//--------------------------------------------------
// This is the ProcessObject that actually does the job:
#include "itkVTKRungeKuttaFinslerIntegration.h"
//--------------------------------------------------
// VTK include to write poly-data:
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
//--------------------------------------------------

#define DIMENSION 3

namespace{
    
    int DoIt( int argc, char * argv[] )
    {
        PARSE_ARGS;
        
        // Only float type is allowed. We DO need such precision:
        typedef float PixelType;
        
        /********************************************************************************************/
        // TYPES FOR THE INPUT/OUTPUT DATA
        // These are the types of images we actually need:
        typedef itk::Image<PixelType,DIMENSION>                                 ScalarImageType;
        typedef itk::Image<unsigned short,DIMENSION>                            LabelImageType;
        typedef itk::CovariantVector<PixelType,DIMENSION>                       CovariantVectorType;
        typedef itk::Image<CovariantVectorType,DIMENSION>                       VectorImageType;
        /********************************************************************************************/
        // TYPES FOR READERS
        // Readers:
        typedef itk::ImageFileReader<ScalarImageType>                           ScalarReaderType;
        typedef ScalarReaderType::Pointer                                       ScalarReaderPointer;
        typedef itk::ImageFileWriter<ScalarImageType>                           ScalarWriterType;
        typedef ScalarWriterType::Pointer                                       ScalarWriterPointer;
        typedef itk::ImageFileReader<LabelImageType>                            LabelReaderType;
        typedef LabelReaderType::Pointer                                        LabelReaderPointer;
        typedef itk::ImageFileReader<VectorImageType>                           VectorReaderType;
        typedef VectorReaderType::Pointer                                       VectorReaderPointer;
        /********************************************************************************************/
        // TYPE FOR THE FILTER
        typedef itk::VTKRungeKuttaFinslerIntegration< VectorImageType,
        ScalarImageType, LabelImageType >   IntegratorType;
        typedef IntegratorType::Pointer                                         IntegratorPointer;
        /********************************************************************************************/
        // TRY TO READ THE IMAGES
        ScalarReaderPointer costsReader      = ScalarReaderType::New();
        LabelReaderPointer  labelReader      = LabelReaderType::New();
        VectorReaderPointer directionsReader = VectorReaderType::New();
        costsReader->SetFileName( inputCost.c_str() );
        labelReader->SetFileName( inputSeeds.c_str() );
        directionsReader->SetFileName( inputDirections.c_str() );
        try{
            costsReader->Update();
            labelReader->Update();
            directionsReader->Update();
        }
        catch( itk::ExceptionObject& e ){
            std::cerr << "I wasn't able to read one or more of the required volumes: "
            << std::endl << e.GetDescription() << std::endl;
            return EXIT_FAILURE;
        }
        /********************************************************************************************/
        // DO FIBER TRACKING:
        IntegratorPointer integrator = IntegratorType::New();
        
        integrator->SetArrivals( directionsReader->GetOutput() );
        integrator->SetCosts( costsReader->GetOutput() );
        integrator->SetLabels( labelReader->GetOutput() );
        
        integrator->SetMinimumLength( iMinimumLength );
        integrator->SetMaximumLength( iMaximumLength );
        integrator->SetStoppingCurvature( iStoppingCurvature );
        integrator->SetStepLength( iStepLength );
        integrator->SetLabel( iLabel );
        
        //------------------
        // Check if we have to output the costs at the target points normalized by
        // the arc length of the computed fiber bundles:
        ScalarImageType::Pointer ncosts;
        if( normalizedCosts.length() > 0 ){
            // We have to. We start by creating a volume with the appropriate size:
            ncosts = ScalarImageType::New();
            ncosts->SetRegions( costsReader->GetOutput()->GetLargestPossibleRegion() );
            ncosts->SetOrigin( costsReader->GetOutput()->GetOrigin() );
            ncosts->SetSpacing( costsReader->GetOutput()->GetSpacing() );
            ncosts->SetDirection( costsReader->GetOutput()->GetDirection() );
            ncosts->Allocate();
            ncosts->FillBuffer( (-2) * itk::NumericTraits<PixelType>::One );
            integrator->SetNormalizedCosts( ncosts );
        }
        else{
            std::cout << "The normalized costs will not be written" << std::endl;
        }
        //------------------
        
        try{
            integrator->Update();
        }
        catch(itk::ExceptionObject& e ){
            std::cerr << "Could not integrate the streamlines: " << e.GetDescription() << std::endl;
            return EXIT_FAILURE;
        }
        /********************************************************************************************/
        // TRY TO WRITE OUT THE OUTPUT (The file has .vtp or .vtk extension):
        vtkSmartPointer<vtkXMLPolyDataWriter> writerVtp = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        vtkSmartPointer<vtkPolyDataWriter>    writerVtk = vtkSmartPointer<vtkPolyDataWriter>::New();
        writerVtp->SetFileName( outputFibers.c_str() );
        writerVtp->SetInput( integrator->GetStreamlines() );
        writerVtk->SetFileName( outputFibers.c_str() );
        writerVtk->SetInput( integrator->GetStreamlines() );
        try{
            if( outputFibers.rfind(".vtp")==outputFibers.size()-4 ){
                if( writerVtp->Write() != 1 ){
                    std::cerr << "I wasn't able to write the output fibers." << std::endl;
                    return EXIT_FAILURE;
                }
            }
            else if( outputFibers.rfind(".vtk")==outputFibers.size()-4 ){
                if( writerVtk->Write() != 1 ){
                    std::cerr << "I wasn't able to write the output fibers." << std::endl;
                    return EXIT_FAILURE;
                }
            }
            else{
                std::cerr << "WARNING: I cannot recognize the file format. Using vtk writer..." << std::endl;
                if( writerVtk->Write() != 1 ){
                    std::cerr << "I wasn't able to write the output fibers." << std::endl;
                    return EXIT_FAILURE;
                }
            }
        }
        catch( ... ){
            std::cerr << "I wasn't able to write the output fibers." << std::endl;
            return EXIT_FAILURE;
        }
        /********************************************************************************************/
        // Try to write out the normalized costs in case it is necessary:
        if( normalizedCosts.length() > 0 ){
            ScalarWriterPointer scwriter = ScalarWriterType::New();
            scwriter->SetInput( integrator->GetNormalizedCosts() );
            scwriter->SetFileName( normalizedCosts.c_str() );
            try{
                scwriter->Update();
            }
            catch( itk::ExceptionObject &e ){
                std::cerr << "Couldn't write to " << normalizedCosts << ": " << std::endl;
                std::cerr << e.GetDescription() << std::endl;
            }
        }
        /********************************************************************************************/
        return EXIT_SUCCESS;
    }
    
}


int main( int argc, char * argv[] )
{
    PARSE_ARGS;
    return DoIt( argc, argv );
}


