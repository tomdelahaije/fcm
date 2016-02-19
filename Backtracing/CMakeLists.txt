cmake_minimum_required(VERSION 2.8)
project( FinslerBacktracing )

################################################################################
set(STANDALONE_CMAKE_HELPER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BuildScripts)
include(${STANDALONE_CMAKE_HELPER_DIR}/CMakeBuildMacros.cmake)
include(${STANDALONE_CMAKE_HELPER_DIR}/PreventInSourceBuilds.cmake)
include(${STANDALONE_CMAKE_HELPER_DIR}/CMakeSTANDALONE3BuildMacros.cmake)

###
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
###
CHECKIFSLICER3BUILD()  ## Call the convenience macro

if (Slicer_SOURCE_DIR)
    if(Slicer_BINARY_DIR)
    else(Slicer_BINARY_DIR)
        # Slicer_BINARY_DIR only exists for Slicer 4. This flag is only used for Slicer3 builds:
        SET_SOURCE_FILES_PROPERTIES(${PROJECT_NAME}.cxx COMPILE_FLAGS -DTHISISASLICERBUILD=1)
    endif(Slicer_BINARY_DIR)
endif (Slicer_SOURCE_DIR)
################################################################################

if(NOT ITK_FOUND)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

if(NOT VTK_FOUND)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif(NOT VTK_FOUND)

###

set(EXTRALIBS
    vtkIO
)


# For Slicer 4 builds, simply call the standard macro:
if(Slicer_BINARY_DIR) # This variable exits only for Slicer4

    slicerMacroBuildCLI(
       NAME ${PROJECT_NAME}
       LOGO_HEADER ${PROJECT_SOURCE_DIR}/ModuleLogo.h
       TARGET_LIBRARIES ${ITK_LIBRARIES} ${EXTRALIBS} ModuleDescriptionParser
       LINK_DIRECTORIES ${ModuleDescriptionParser_BINARY_DIR}
       INCLUDE_DIRECTORIES ${SlicerBaseCLI_SOURCE_DIR} ${Slicer_SOURCE_DIR}/Applications/CLI/DiffusionApplications/DiffusionApplicationsCommon
       ${SlicerBaseCLI_BINARY_DIR}
    )
    # This is necessary to tell the main .cxx what kind of build we are using:
    SET_SOURCE_FILES_PROPERTIES(${PROJECT_NAME}.cxx COMPILE_FLAGS -DSLICERV4=1)
    
else(Slicer_BINARY_DIR)
 
    set(CLP FinslerBacktracing)
    set ( ${CLP}_SOURCE ${CLP}.cxx )
    CONFIGURESTANDALONEORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "${EXTRALIBS}" "")
    target_link_libraries( ${CLP} ITKIO ITKAlgorithms )
    
    INSTALL(TARGETS FinslerBacktracing DESTINATION bin)
    
endif(Slicer_BINARY_DIR)

# The below is for testing purposes. The first 6 tests (i.e. the first three couples of tests) are for checking the parallel fast-sweeping code. The two additional tests are for checking the HARDI code


#Testing:
if (BUILD_TESTING)

    add_executable(${CLP}_Compare ${CLP}_Compare.cxx)
    target_link_libraries(${CLP}_Compare
        vtkIO
    )
    
    set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/../TestData)
    if (Slicer_SOURCE_DIR)
        
        # This is a Slicer build
        # Depending on the actual version (3 or 4), the calling sequence differs:
        if(Slicer3_BINARY_DIR)
            # This is for Slicer3
            set(OUTPUT_DIR ${Slicer3_BINARY_DIR}/Testing)
            set(Slicer_EXE ${Slicer3_BINARY_DIR}/Slicer3)
            SET(RUNEXECB
                ${Slicer_EXE}
                --launch
                ${CLP}
            )
            SET(COMPEXECB
                ${Slicer_EXE}
                --launch
                ${CLP}_Compare
            )
        endif(Slicer3_BINARY_DIR)
        
        if(Slicer_BINARY_DIR)
            # And this is for Slicer4. It requires calling the CLI with the full path
            set(OUTPUT_DIR ${Slicer_BINARY_DIR}/Testing)
            set(Slicer_EXE ${Slicer_BINARY_DIR}/Slicer)
            SET(RUNEXECB
                ${Slicer_EXE}
                --launch
                "${Slicer_BINARY_DIR}/lib/Slicer-4.1/cli-modules/${CLP}"
            )
            SET(COMPEXECB
                ${Slicer_EXE}
                --launch
                "${Slicer_BINARY_DIR}/bin/${CLP}_Compare"
            )
        endif(Slicer_BINARY_DIR)
        
    else( Slicer_SOURCE_DIR )
        set(OUTPUT_DIR ${CMAKE_BINARY_DIR}/Testing)
        SET( RUNEXECB ${CMAKE_BINARY_DIR}/${CLP} )
        SET( COMPEXECB ${CMAKE_BINARY_DIR}/${CLP}_Compare )
    endif (Slicer_SOURCE_DIR)
    
    add_test(${PROJECT_NAME}_SampleTest
        ${RUNEXECB}
        --minimumlength 10
        --maximumlength 800
        --stoppingCurvature 0.8
        --stepLength 0.5
        --label 1
        ${OUTPUT_DIR}/sampleDirections.nrrd
        ${TEST_DATA_DIR}/sampleFinslerCosts.nii.gz
        ${TEST_DATA_DIR}/sampleSeeds.nii.gz
        ${OUTPUT_DIR}/sampleFibers.vtp
    )
    
    add_test(${PROJECT_NAME}_SampleCompare
        ${COMPEXECB}
        ${TEST_DATA_DIR}/sampleFibers.vtp
        ${OUTPUT_DIR}/sampleFibers.vtp
    )
    
    add_test(${PROJECT_NAME}_CleanResults
        rm
        ${OUTPUT_DIR}/sampleDirections.nrrd
        ${OUTPUT_DIR}/sampleFibers.vtp
    )
  
endif (BUILD_TESTING)
