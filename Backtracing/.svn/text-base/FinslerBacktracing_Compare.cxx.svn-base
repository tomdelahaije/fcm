/*=========================================================================
 
 Program:   Diffusion Applications
 Language:  C++
 Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/dwiNoiseFilter/dwiNoiseFilter.cxx $
 Date:      $Date: 2008-11-25 18:46:58 +0100 (Tue, 25 Nov 2008) $
 Version:   $Revision: 7972 $
 
 Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.
 
 See License.txt or http://www.slicer.org/copyright/copyright.txt for details.
 
 ==========================================================================*/

// This is a dumb program to compare two VTK polydata objects
// stored in two files, returning success if they are equal
// and failure otherwise
#include <math.h>
#include "vtkXMLPolyDataReader.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"



#define DIMENSION 3



int main( int argc, const char * argv[] )
{
    if( argc<3 ){
        std::cerr << "Usage: " << argv[0] << " xmlPolyData1.vtp xmlPolyData2.vtp" << std::endl;
        return EXIT_FAILURE;
    }
    
    vtkSmartPointer<vtkXMLPolyDataReader> reader1 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    vtkSmartPointer<vtkXMLPolyDataReader> reader2 = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    
    reader1->SetFileName( argv[1] );
    reader2->SetFileName( argv[2] );
    
    if( !reader1->CanReadFile( argv[1] ) ){
        std::cerr << "Unable to read: " << argv[1] << " with vtkXMLPolyDataReader" << std::endl;
        return EXIT_FAILURE;
    }
    
    if( !reader2->CanReadFile( argv[2] ) ){
        std::cerr << "Unable to read: " << argv[2] << " with vtkXMLPolyDataReader" << std::endl;
        return EXIT_FAILURE;
    }
    
    try
    {
        reader1->Update();
    }
    catch ( ... )
    {
        std::cerr << "Unable to read: " << argv[1] << " with vtkXMLPolyDataReader" << std::endl;
        return EXIT_FAILURE;
    }
    
    try
    {
        reader2->Update();
    }
    catch ( ... )
    {
        std::cerr << "Unable to read: " << argv[2] << " with vtkXMLPolyDataReader" << std::endl;
        return EXIT_FAILURE;
    }
    
    vtkSmartPointer<vtkPolyData> polydata1 = reader1->GetOutput();
    vtkSmartPointer<vtkPolyData> polydata2 = reader2->GetOutput();
    
    if( polydata1==NULL ){
        std::cerr << "NULL pinter from reader1!" << std::endl;
        return  EXIT_FAILURE;
    }
    
    if( polydata2==NULL ){
        std::cerr << "NULL pinter from reader2!" << std::endl;
        return  EXIT_FAILURE;
    }
    
    // Now we have the polt data. We can start comparing both fiber-bundles:
    if( polydata1->GetNumberOfPoints() == polydata2->GetNumberOfPoints() ){
        std::cout << "There are " << polydata1->GetNumberOfPoints() << " points in each model" << std::endl;
        if( polydata1->GetNumberOfCells() == polydata2->GetNumberOfCells() ){
            std::cout << "There are " << polydata1->GetNumberOfCells() << " streamlines in each model" << std::endl;
            for( vtkIdType cc=0; cc<polydata1->GetNumberOfCells(); ++cc ){
                vtkSmartPointer<vtkIdList> cell1 = vtkSmartPointer<vtkIdList>::New();
                vtkSmartPointer<vtkIdList> cell2 = vtkSmartPointer<vtkIdList>::New();
                polydata1->GetCellPoints( cc, (vtkIdList*)cell1 );
                polydata2->GetCellPoints( cc, (vtkIdList*)cell2 );
                if( cell1->GetNumberOfIds() == cell2->GetNumberOfIds() ){
                    for( vtkIdType ii=0; ii<cell1->GetNumberOfIds(); ++ii ){
                        if( cell1->GetId(ii) == cell2->GetId(ii) ){
                            double p1[3];
                            double p2[3];
                            polydata1->GetPoint( cell1->GetId(ii), (double*)p1 );
                            polydata2->GetPoint( cell2->GetId(ii), (double*)p2 );
                            double dist  = ( p1[0]-p2[0] )*( p1[0]-p2[0] );
                            dist        += ( p1[1]-p2[1] )*( p1[1]-p2[1] );
                            dist        += ( p1[2]-p2[2] )*( p1[2]-p2[2] );
                            if( dist>1e-3 ){
                                std::cerr << "Different " << ii << "-th point in streamline " << cc << std::endl;
                                std::cerr << "Distance is: " << dist << std::endl;
                                return EXIT_FAILURE;
                            }
                        }
                        else{
                            std::cerr << "Different " << ii << "-th identifier in streamline " << cc << std::endl;
                            return EXIT_FAILURE;
                        }
                    }
                }
                else{
                    std::cerr << "Different number of points in cell: " << cc << std::endl;
                    std::cerr << "      model1: " << cell1->GetNumberOfIds() << std::endl;
                    std::cerr << "      model2: " << cell2->GetNumberOfIds() << std::endl;
                    return EXIT_FAILURE;
                }
            }
        }
        else{
            std::cerr << "Different number of streamlines in each polyData" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else{
        std::cerr << "Different number of points in each polyData" << std::endl;
        return EXIT_FAILURE;
    }
        
    return EXIT_SUCCESS;
}


