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

#ifdef THISISASLICERBUILD

#include "itkPluginUtilities.h"

#else

#ifdef SLICERV4
#include "itkPluginUtilities.h"
#else
#include "SlicerExecutionModel/itkPluginUtilities.h"
#endif

#endif

#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include "itkCastImageFilter.h"
#include "FinslerTractographyCLP.h"
//--------------------------------------------------
// Specific includes:
#include "itkParallelFastSweeping.h"
#include "HARDI-ITK/sphericalHarmonics.h"
#include "HARDI-ITK/itkComputeSHCoefficientsFilter.h"
#include "HARDI-ITK/itkComputeLocalCostFromSH.h"
#include "itkChangeLabelImageFilter.h"
#include "itkExpNegativeImageFilter.h"
#include "Generate60Gradients.h"
//--------------------------------------------------
#define DIMENSION 3


namespace{
    
    template<class PixelType>
    int DoIt( int argc, char * argv[], PixelType )
    {
        PARSE_ARGS;
        // do the typedefs
        typedef itk::VectorImage<PixelType,DIMENSION>                            DiffusionImageType;
        typedef itk::Image<float,DIMENSION>                                      ScalarImageType;
        typedef itk::Image<unsigned short,DIMENSION>                             LabelImageType;
        typedef itk::Image<itk::CovariantVector<float,DIMENSION>,DIMENSION>      VectorImageType;
        
        typedef itk::CovariantVector<double,DIMENSION>                           CovariantVectorType;
        
        typedef itk::ImageFileReader<LabelImageType>                             LabelReaderType;
        
        typedef itk::VectorImage<float,DIMENSION>                                SHImageType;      // Each voxel is the list of SH coefficients
        typedef itk::VectorImage<float,DIMENSION>                                CostsImageType;   // Each voxel is the directional cost
        // Filter to compute the SH coefficients from the DWI:
        typedef itk::ComputeSHCoefficientsFilter<DiffusionImageType,SHImageType> SHComputeType;
        typedef typename SHComputeType::Pointer                                  SHComputePointer;
        typedef typename SHComputeType::GradientType                             GradientType;
        typedef typename SHComputeType::ListOfGradientsType                      ListOfGradientsType;
        // Filter to compute the local costs from the SH coefficients
        typedef itk::ComputeLocalCostFromSH<SHImageType,CostsImageType>          CostsComputeType;
        typedef typename CostsComputeType::Pointer                               CostsComputePointer;
        
        typedef itk::ParallelFastSweeping<CostsImageType,ScalarImageType,LabelImageType> FinslerDistanceComputeType;
        typedef typename FinslerDistanceComputeType::Pointer                             FinslerDistanceComputePointer;
        
        // This filter is used to remove Inf pixels that arise when the mask has white islands
        // without seeding points. Such pixels do not affect the way Fast Sweeping works (this
        // scenario is accounted for in the program logic), but can be quite anoying when we
        // view the costs map because the window/level controls go crazy. We simply substitute
        // Infinite values with negative values:
        typedef itk::ChangeLabelImageFilter<ScalarImageType,ScalarImageType>     RelabelFilterType;
        typedef RelabelFilterType::Pointer                                       RelabelFilterPointer;
        // This filter is used only if the local costs are computed as minus de logarithm
        // of the OPDF:
        typedef itk::ExpNegativeImageFilter<ScalarImageType,ScalarImageType>     ExponentialFilterType;
        typedef ExponentialFilterType::Pointer                                   ExponentialFilterPointer;
        
        typedef itk::ImageFileWriter<ScalarImageType> ScalarWriterType;
        typedef itk::ImageFileWriter<VectorImageType> VectorWriterType;
        
        //================================================================================
        //================================================================================
        /** I- READ THE INPUT DWI VOLUME AND INTERPRET THE METADATA */
        
        std::vector< CovariantVectorType >                   diffusionDirections;
        
        typedef itk::ImageFileReader<DiffusionImageType>     FileReaderType;
        typename FileReaderType::Pointer reader = FileReaderType::New();
        reader->SetFileName( inputDWI.c_str() );
        try{
            reader->Update();
        }
        catch( itk::ExceptionObject& e ){
            std::cerr << "Cannot read DWI volume:" << std::endl << e.GetDescription() << std::endl;
            return EXIT_FAILURE;
        }
        
        typedef itk::MetaDataDictionary DictionaryType;
        const DictionaryType & dictionary = reader->GetMetaDataDictionary();
        
        typedef itk::MetaDataObject< std::string >                       MetaDataStringType;
        typedef itk::MetaDataObject< std::vector<std::vector<double> > > MetaDataMFrameType;
        DictionaryType::ConstIterator itr = dictionary.Begin();
        DictionaryType::ConstIterator end = dictionary.End();
        
        double       dBValue                = 1000;
        bool         iFoundBValue           = false;
        unsigned int channels               = 0;
        //%---------------------------------------------------------------------------
        // We also need to parse the space directions and the measurement frame; this is
        // because we need to relate the gradient directions to directions in the
        // ijk-space to appropriately compute the costs.
        //
        //   FOR IMPORTANT INFORMATION ON THIS, LOOK:
        //   http://www.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
        //
        /* [space directions]: The "space directions" field defines, one column vector at a
         time, the transform from IJK raster image coordinates (I==fastest;K==slowest) to RAS
         anatomical image coordinates.
         */
        /* [measurement frame]: This defines, one column vector at a time, the matrix which
         transforms the (X,Y,Z) coordinates used for listing the diffusion-sensitizing
         gradient directions (below), into whatever space is named by "space", in this case
         RAS anatomical image coordinates.
         */
        // It's irrelevant if the space is RAS or LPS, since the measurement frame is
        // defined regardless of the actual space.
        //%---------------------------------------------------------------------------
        typename DiffusionImageType::DirectionType measurementFrame;
        measurementFrame.SetIdentity();
        typename DiffusionImageType::DirectionType spaceDirections = reader->GetOutput()->GetDirection();
        bool     iFoundMeasurementFrame = false;
        
        while( itr != end ){
            // Make sure this entry of the dictionary is not empty:
            itk::MetaDataObjectBase::Pointer entry = itr->second;
            if( entry.GetPointer() ){
                // Examine the tag key to check if it is one of those we are interested in
                ::size_t pos = itr->first.find("DWMRI_gradient"); // Is it a gradient direction?
                if ( pos != std::string::npos ){
                    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>( entry.GetPointer() );
                    std::string tagkey   = itr->first;
                    std::string tagvalue = entryvalue->GetMetaDataObjectValue();
                    double dx[DIMENSION];
                    std::sscanf(tagvalue.c_str(), "%lf %lf %lf\n", &dx[0], &dx[1], &dx[2]);
                    CovariantVectorType dGrad = (CovariantVectorType)(dx);
                    /** We need to be able to cope with nrrd files for which different
                     b-values are encoded as the norm of the diffusion gradient. Hence, we
                     remove the normalization:
                     if( dGrad.GetNorm()>1e-6 )
                     dGrad.Normalize();
                     */
                    diffusionDirections.push_back( dGrad );
                    ++channels;
                }
                // Examine the tag key to check if it is one of those we are interested in
                pos = itr->first.find("DWMRI_b-value"); // Is it the b-value?
                if ( pos != std::string::npos ){
                    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>( entry.GetPointer() );
                    std::string tagvalue = entryvalue->GetMetaDataObjectValue();
                    std::sscanf(tagvalue.c_str(), "%lf\n", &dBValue );
                    iFoundBValue = true;
                }
                // Examine the tag key to check if it is one of those we are interested in
                pos = itr->first.find("NRRD_measurement frame"); // Is it the measurement frame?
                if ( pos != std::string::npos ){
                    MetaDataMFrameType::Pointer entryvalue = dynamic_cast<MetaDataMFrameType *>( entry.GetPointer() );
                    std::vector< std::vector<double> > tagvalue = entryvalue->GetMetaDataObjectValue();
                    // The measurement frame is stored in the tagvalue one column at a time, i.e., in
                    // a nrrd file we have [[m00,m10,m20],[m01,m11,m21],[m02,m12,m22]]
                    for( unsigned int cm=0; cm<3; ++cm ){
                        // This is each COLUMN of the measurement frame
                        for( unsigned int rm=0; rm<3; ++rm ){
                            // This is each ROW of the measurement frame
                            measurementFrame[rm][cm] = tagvalue[cm][rm];
                        }
                    }
                    iFoundMeasurementFrame = true;
                }
            } // END if( entryvalue )
            ++itr;
        } // END while( itr != end )
        
        if( !iFoundBValue )
            std::cerr << "WARNING: I didn't find a b-value. Using: " << dBValue << std::endl;
        if( !iFoundMeasurementFrame )
            std::cerr << "WARNING: I didn't find a measurement frame. Using identity" << std::endl;
        //%---------------------------------------------------------------------------
        
        // Find the first zero baseline image and use it for the noise estimation
        ::size_t iNrOfDWIs      = diffusionDirections.size();
        ::size_t iFirstBaseline = std::string::npos;
        for ( ::size_t iI=0; iI<iNrOfDWIs; iI++ ){
            if ( diffusionDirections[iI].GetNorm() < 1e-6 ){
                iFirstBaseline = iI;
                break;
            }
        }
        if ( iFirstBaseline == std::string::npos ){
            std::cout << "Did not find an explicit baseline image." << std::endl;
            std::cout << "Treating the first volume as the baseline volume." << std::endl;
            iFirstBaseline = 0;
        }
        //================================================================================
        //================================================================================
        /** II- READ THE VOLUME OF SEEDING POINTS */
        typename LabelReaderType::Pointer seedsReader = LabelReaderType::New();
        seedsReader->SetFileName( inputSeeds.c_str() );
        try{
            seedsReader->Update();
        }
        catch( itk::ExceptionObject& e ){
            std::cerr << "Cannot read seeds volume:" << std::endl << e.GetDescription() << std::endl;
            return EXIT_FAILURE;
        }
        //================================================================================
        //================================================================================
        /** III- READ THE VOLUME WITH THE MASK */
        typename LabelReaderType::Pointer maskReader = LabelReaderType::New();
        maskReader->SetFileName( inputMask.c_str() );
        bool useMask;
        if( inputMask.length() > 0 ){
            try{
                maskReader->Update();
                // Check if this has the appropriate size
                if( maskReader->GetOutput()->GetLargestPossibleRegion() != reader->GetOutput()->GetLargestPossibleRegion() ){
                    std::cout << "The mask has a weird size. I'm not using it..." << std::endl;
                    useMask = false;
                }
                else
                    useMask = true;
            }
            catch( itk::ExceptionObject& e ){
                std::cout << "I dindn't find a mask. Running without it..." << std::endl;
                useMask = false;
            }
        }
        else{
            std::cout << "I'm not using the mask" << std::endl;
            useMask = false;
        }
        //================================================================================
        //================================================================================
        /** IV- COMPUTE THE VOLUME OF SH AND THE DIRECTIONAL COSTS */
        // Create filters and do the job:
        SHComputePointer              shCompute      = SHComputeType::New();
        CostsComputePointer           costsCompute   = CostsComputeType::New();
        // SET THE PARAMETERS AND THE INPUT TO THE FILTER COMPUTING THE SH COEFFICIENTS
        // Set the list of gradient directions; for the algorithm to work, these
        // gradients must refer to the ijk space, so we need to convert those
        // retrieved from the DWI volume to the correct space:
        std::vector< CovariantVectorType > ijkDirs( diffusionDirections.size() );
        typename DiffusionImageType::DirectionType grad2ijk;
        grad2ijk  = spaceDirections.GetInverse();
        grad2ijk *= measurementFrame;
        // This vector is used to store the b-value for each channel:
        std::vector<double> dBValues_list( diffusionDirections.size() );
        for( unsigned int l=0; l<ijkDirs.size(); ++l ){
            // ----- Old:
            // CovariantVectorType correctedGrad = grad2ijk * diffusionDirections[l];
            // ----- New:
            CovariantVectorType correctedGrad = diffusionDirections[l];
            // The gradients can be unnormalized, standing for a different b-value in
            // each channel. Hence, we need to normalize the gradients and keep track
            // of the corresponding b-values:
            dBValues_list[l] = ( correctedGrad.GetNorm() * dBValue );
            correctedGrad.Normalize();
            correctedGrad = grad2ijk * correctedGrad;
            // ----- END new
            correctedGrad.Normalize();
            ijkDirs[l] = correctedGrad;
        }
        shCompute->SetList( ijkDirs );
        // Set the b-value of the acquisition (retireve from the DWI volume):
        shCompute->SetBValues( dBValues_list );
        // Set the "lambda" regularization parameter (retrieve from the GUI):
        shCompute->SetLambda( iLambda );
        // Set the degree of the SH expansions (retireve from the GUI;
        //       future work: check or automatically determine from the number of DWI channels):
        shCompute->SetL( iLSH );
        // Set the DWI-related measure according to the local cost to be computed:
        if( iTypeOfLocalCost.compare("1/ADC/Psi^3") == 0 )
            shCompute->SetComputeATTforOPDT();
        else if( iTypeOfLocalCost.compare("E(t|R)") == 0 ){
            shCompute->SetComputeMeanTime();
            std::cout << "I'm using the experimental local directional cost based on E(t|R)" << std::endl;
        }
        else
            shCompute->SetComputeATT();
        // Set the input:
        shCompute->SetInput( reader->GetOutput() );
        //-----------------------------------------------------------------------------------------------------------------------------------
        // SET THE PARAMETERS AND THE INPUT TO THE FILTER COMPUTING THE LOCAL COST
        // Default neighboring directions (Future work: accept these directions as a parameter):
        typedef FinslerDistanceComputeType::DirectionType NeighborDirectionType;
        std::vector<NeighborDirectionType> ndirs;
        NeighborDirectionType              ndir;
        if( iNumDirs.compare("26") == 0 ){
            for( int z=-1; z<=1; ++z ){
                for( int y=-1; y<=1; ++y ){
                    for( int x=-1; x<=1; ++x ){
                        if( x!=0 || y!=0 || z!=0 ){
                            ndir[0] = x;
                            ndir[1] = y;
                            ndir[2] = z;
                            ndir.Normalize();
                            ndirs.push_back( ndir );
                        }
                    }
                }
            }
        }
        else if( iNumDirs.compare("60") == 0 ){
            const unsigned int rows = 30;
            const unsigned int cols = 3;
            ndirs.resize(rows);
            Generate60Gradients< std::vector<NeighborDirectionType> >( ndirs, rows, cols );
            for( unsigned int k=0; k<rows; ++k ){
                ndirs[k].Normalize();
                ndirs.push_back( -ndirs[k] );
            }
        }
        else if( iNumDirs.compare("all") == 0 ){
            for( unsigned int k=0; k<diffusionDirections.size(); ++k ){
                CovariantVectorType dGrad = diffusionDirections[k];
                if( dGrad.GetNorm()>1e-6 ){
                    dGrad.Normalize(); // Before this point, gradients are not normalized
                    ndir[0] = dGrad[0];
                    ndir[1] = dGrad[1];
                    ndir[2] = dGrad[2];
                    ndirs.push_back(  ndir );
                    ndirs.push_back( -ndir );
                }
            }
        }
        else{
            std::cerr << "Unrecognized option: " << iNumDirs << std::endl;
            return EXIT_FAILURE;
        }
        
        // Set the neighboring directions:
        costsCompute->SetNeighboringDirections( ndirs );
        // How should we compute the local costs?
        if( iTypeOfLocalCost.compare("(E(q)/Phi(r))^3") == 0 )
            costsCompute->SetComputeATTOPDF();
        else if( iTypeOfLocalCost.compare("1/Psi(r)^3") == 0 )
            costsCompute->SetComputeOPDFINV();
        else if( iTypeOfLocalCost.compare("1/ADC/Psi^3") == 0 ){
            costsCompute->SetComputeDIFFUSIVITY();
        }
        else if( iTypeOfLocalCost.compare("E(t|R)") == 0 ){
            costsCompute->SetComputeMEANTIME();
            std::cout << "I'm using the experimental local directional cost based on E(t|R)" << std::endl;
        }
        else if( iTypeOfLocalCost.compare("-log(Psi(r))") == 0 )
            costsCompute->SetComputeOPDFLOG();
        else if( iTypeOfLocalCost.compare("1") == 0 ){
            costsCompute->SetComputeISOTROPIC();
            std::cout << "I'm using constant and isotropic local directional cost" << std::endl;
        }
        else if( iTypeOfLocalCost.compare("test") == 0 )
            costsCompute->SetComputeTEST();
        else{
            std::cerr << "Unsupported type of local cost: " << iTypeOfLocalCost << std::endl;
            return EXIT_FAILURE;
        }
        // Set the input:
        costsCompute->SetInput( shCompute->GetOutput() );
        //================================================================================
        //================================================================================
        /** V- DO FAST SWEEPING */
        FinslerDistanceComputePointer fastSweeping = FinslerDistanceComputeType::New();
        fastSweeping->SetInput( costsCompute->GetOutput() );
        //-----------------------------------------------------------------------------------------------------------------------------------
        // SET THE PARAMETERS AND THE INPUT TO THE FILTER PERFORMING FAST SWEEPING
        // Set the maximum number of iterations (retrieve from the GUI):
        fastSweeping->SetMaxIters( iMaxIters );
        // Set the cost factor for the stop criterion (retrieve log-value from the GUI):
        double factor = (double)iCostFraction;
        factor = ::exp( factor * ::log(10.0f) );
        fastSweeping->SetCostFactor( factor );
        // Set the neighboring directions:
        fastSweeping->SetNeighboringDirections( ndirs );
        // Set the acceleration parametes:
        fastSweeping->SetUseAcceleration( iAcceleration );
        fastSweeping->SetAccelerateIter( iAccelerateIter );
        // Check if multi-threading is to be used:
        // Update: Nov. 2014: threading disabled until threading functionality is fixed
        fastSweeping->SetUseThreads( 0 /* iUseThreads */ );
        // Set the input
        fastSweeping->SetInput( costsCompute->GetOutput() );
        //Set the output file name
        fastSweeping->SetOutputConnMapName( outputConnMap.c_str() );
        // Set the seeding points:
        std::vector<ScalarImageType::IndexType > seeds;
        seeds.clear();
        typedef itk::ImageRegionConstIteratorWithIndex<LabelImageType> IteratorType;
        IteratorType it( seedsReader->GetOutput(), seedsReader->GetOutput()->GetLargestPossibleRegion() );
        for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
            if( it.Get() == iLabel ){
                // This is an index in the mask image with the seeing points. It does not
                // necessarily correspond to the same index in the ODF/DWI image, so we have
                // to double convert from index to physical coord in the mask and from physical
                // coord to index in the ODF image:
                typename DiffusionImageType::PointType physical;
                typename DiffusionImageType::IndexType indexed = it.GetIndex();
                // Transform the index in the labelmap to a physical point:
                seedsReader->GetOutput()->TransformIndexToPhysicalPoint( indexed, physical );
                // Transform the physical point obtained to an index in the DWI volume:
                bool amiinside = reader->GetOutput()->TransformPhysicalPointToIndex	( physical, indexed );
                // Of the point is inside the DWI volume buffer, set the new seeding point:
                if( amiinside )
                    seeds.push_back( indexed );
                else
                    std::cout << "WARNING: The index " << it.GetIndex() << " maps to an index " << indexed << " outside the image" << std::endl; 
            }
        }
        // Make sure at least one seed is present:
        if( seeds.size()==0 ){
            std::cerr << "I have no seeds to work with!" << std::endl;
            return EXIT_FAILURE;
        }
        // Set the seeding points of the Finsler filter:
        fastSweeping->SetSeedingPoints( seeds );
        //-----------------------------------------------------------------------------------------------------------------------------------
        // Set the mask
        if( useMask )
            fastSweeping->SetMask( maskReader->GetOutput() );
        //-----------------------------------------------------------------------------------------------------------------------------------
        // Avoid infinites
        RelabelFilterPointer relabelFilter = RelabelFilterType::New();
        relabelFilter->SetInput( fastSweeping->GetOutput() );
        relabelFilter->SetChange( itk::NumericTraits<ScalarImageType::PixelType>::max(), static_cast<ScalarImageType::PixelType>(-1) );
        //-----------------------------------------------------------------------------------------------------------------------------------
        // Prepare to write
        typename ScalarWriterType::Pointer writer = ScalarWriterType::New();
        writer->SetFileName( outputCost.c_str() );
        //-----------------------------------------------------------------------------------------------------------------------------------
        // Compute the exponential if needed
        ExponentialFilterPointer exponentialFilter;
        if( iTypeOfLocalCost.compare("-log(Psi(r))") == 0 ){
            exponentialFilter = ExponentialFilterType::New();
            exponentialFilter->SetInput( relabelFilter->GetOutput() );
            exponentialFilter->SetFactor( 1.0f );
            // Taking the exponential does not seem to work well. By now, we keep the log-value:
            writer->SetInput( relabelFilter->GetOutput() );
        }
        else
            writer->SetInput( relabelFilter->GetOutput() );
        //-----------------------------------------------------------------------------------------------------------------------------------
        // Update the whole pipeline:
        writer->UseInputMetaDataDictionaryOff();
        writer->SetMetaDataDictionary( reader->GetMetaDataDictionary() );
        try{
            writer->Update();
        }
        catch( itk::ExceptionObject& e ){
            std::cerr << "Cannot update the filter, something went wrong..." << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            return EXIT_FAILURE;
        }
        //================================================================================
        //================================================================================
        /** VI- WRITE OUT THE MAP OF ARRIVAL DIRECTIONS */
        if( outputDirections.length() > 0 ){
            typename VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
            vectorWriter->SetFileName( outputDirections.c_str() );
            vectorWriter->SetInput( fastSweeping->GetOptimumDirectionsMap() );
            vectorWriter->UseInputMetaDataDictionaryOff();
            vectorWriter->SetMetaDataDictionary( reader->GetMetaDataDictionary() );
            try{
                vectorWriter->Update();
                //std::cout << fastSweeping << std::endl;
            }
            catch( itk::ExceptionObject& e ){
                std::cout << "I cannot write the map of optimal output directions" << std::endl;
            }
        }
        else{
            std::cout << "I am not writting the map of optimum directions" << std::endl;
        }
        //================================================================================
        //================================================================================
        return EXIT_SUCCESS;
    }
    
}

int main( int argc, char * argv[] )
{
    PARSE_ARGS;
    
    itk::ImageIOBase::IOPixelType pixelType;
    itk::ImageIOBase::IOComponentType componentType;
    itk::GetImageType (inputDWI, pixelType, componentType);
    
    // This filter handles all types
    
    switch (componentType)
    {
#ifndef WIN32
        case itk::ImageIOBase::UCHAR:
            return DoIt( argc, argv, static_cast<unsigned char>(0));
            break;
        case itk::ImageIOBase::CHAR:
            return DoIt( argc, argv, static_cast<char>(0));
            break;
#endif
        case itk::ImageIOBase::USHORT:
            return DoIt( argc, argv, static_cast<unsigned short>(0));
            break;
        case itk::ImageIOBase::SHORT:
            return DoIt( argc, argv, static_cast<short>(0));
            break;
        case itk::ImageIOBase::UINT:
            return DoIt( argc, argv, static_cast<unsigned int>(0));
            break;
        case itk::ImageIOBase::INT:
            return DoIt( argc, argv, static_cast<int>(0));
            break;
#ifndef WIN32
        case itk::ImageIOBase::ULONG:
            return DoIt( argc, argv, static_cast<unsigned long>(0));
            break;
        case itk::ImageIOBase::LONG:
            return DoIt( argc, argv, static_cast<long>(0));
            break;
#endif
        case itk::ImageIOBase::FLOAT:
            return DoIt( argc, argv, static_cast<float>(0));
            break;
        case itk::ImageIOBase::DOUBLE:
            std::cout << "DOUBLE type not currently supported." << std::endl;
            return EXIT_FAILURE;
            break;
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
            std::cout << "unknown component type" << std::endl;
            return EXIT_FAILURE;
            break;
    }
    
    return EXIT_SUCCESS;
}


