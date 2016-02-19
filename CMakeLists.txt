cmake_minimum_required(VERSION 2.8)
project( FinslerTractography )

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
else( Slicer_SOURCE_DIR )
    ADD_SUBDIRECTORY(SlicerExecutionModel)
    OPTION(BUILD_TESTING "Build tests for this program" ON)
    if( BUILD_TESTING )
        ENABLE_TESTING()
    endif( BUILD_TESTING )
endif (Slicer_SOURCE_DIR)
################################################################################


OPTION(BUILD_BACKTRACING "Build also the backtracing code (needs VTK)" ON )
if( BUILD_BACKTRACING )
    ADD_SUBDIRECTORY( Backtracing )
endif( BUILD_BACKTRACING )


if(NOT ITK_FOUND)
    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})
endif(NOT ITK_FOUND)

###
set( ADDITIONAL_SOURCES
    HARDI-ITK/sh2hot.cxx
    HARDI-ITK/sphericalHarmonics.cxx
)


# For Slicer 4 builds, simply call the standard macro:
if(Slicer_BINARY_DIR) # This variable exits only for Slicer4

    slicerMacroBuildCLI(
       NAME ${PROJECT_NAME}
       LOGO_HEADER ${PROJECT_SOURCE_DIR}/ModuleLogo.h
       ADDITIONAL_SRCS ${ADDITIONAL_SOURCES}
       TARGET_LIBRARIES ${ITK_LIBRARIES} ModuleDescriptionParser
       LINK_DIRECTORIES ${ModuleDescriptionParser_BINARY_DIR}
       INCLUDE_DIRECTORIES ${SlicerBaseCLI_SOURCE_DIR} ${Slicer_SOURCE_DIR}/Applications/CLI/DiffusionApplications/DiffusionApplicationsCommon
       ${SlicerBaseCLI_BINARY_DIR}
    )
    # This is necessary to tell the main .cxx what kind of build we are using:
    SET_SOURCE_FILES_PROPERTIES(${PROJECT_NAME}.cxx COMPILE_FLAGS -DSLICERV4=1)
    
else(Slicer_BINARY_DIR)
 
    # Include the .txx files in the project so that they are in the xcode, visual studio... workspace
    set(EXTRA_HEADERS
        itkImageDirectionalConstIteratorWithIndex.txx
        itkImageDirectionalIteratorWithIndex.txx
        itkIndexToDirectionImageFilter.txx
        itkParallelFastSweeping.txx
        itkParallelFastSweepingStep.txx
        HARDI-ITK/itkComputeLocalCostFromSH.txx
        HARDI-ITK/itkComputeSHCoefficientsFilter.txx
    )
    set(CLP FinslerTractography)
    set ( ${CLP}_SOURCE ${CLP}.cxx ${ADDITIONAL_SOURCES})
    CONFIGURESTANDALONEORSLICERPROPERTIES(${CLP} ${CLP}.xml "${${CLP}_SOURCE}" "" "${EXTRA_HEADERS}")
    target_link_libraries( ${CLP} ${ITK_LIBRARIES})
    
    INSTALL(TARGETS FinslerTractography DESTINATION bin)
 
endif(Slicer_BINARY_DIR)

# The below is for testing purposes. The first 6 tests (i.e. the first three couples of tests) are for checking the parallel fast-sweeping code. The two additional tests are for checking the HARDI code

if (BUILD_TESTING)
  add_executable(${CLP}_Compare ${CLP}_Compare.cxx)
  target_link_libraries(${CLP}_Compare
    ITKIO
    ITKNumerics
    ITKStatistics
  )
  
  set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/TestData)
  
  if (Slicer_SOURCE_DIR)
  
    # This is a Slicer build
    # Depending on the actual version (3 or 4), the calling sequence differs:
    if(Slicer3_BINARY_DIR)
        # This is for Slicer3
        set(OUTPUT_DIR ${Slicer3_BINARY_DIR}/Testing)
        set(Slicer_EXE ${Slicer3_BINARY_DIR}/Slicer3)
        SET(RUNEXEC
            ${Slicer_EXE}
            --launch
            ${CLP}
        )
        SET(COMPEXEC
            ${Slicer_EXE}
            --launch
            ${CLP}_Compare
        )
    endif(Slicer3_BINARY_DIR)
    
    if(Slicer_BINARY_DIR)
        # And this is for Slicer4. It requires calling the CLI with the full path
        set(OUTPUT_DIR ${Slicer_BINARY_DIR}/Testing)
        set(Slicer_EXE ${Slicer_BINARY_DIR}/Slicer)
        SET(RUNEXEC
            ${Slicer_EXE}
            --launch
            "${Slicer_BINARY_DIR}/lib/Slicer-4.1/cli-modules/${CLP}"
        )
        SET(COMPEXEC
            ${Slicer_EXE}
            --launch
            "${Slicer_BINARY_DIR}/bin/${CLP}_Compare"
        )
    endif(Slicer_BINARY_DIR)
    
  else( Slicer_SOURCE_DIR )
    set(OUTPUT_DIR ${CMAKE_BINARY_DIR}/Testing)
    SET( RUNEXEC ${CMAKE_BINARY_DIR}/${CLP} )
    SET( COMPEXEC ${CMAKE_BINARY_DIR}/${CLP}_Compare )
  endif (Slicer_SOURCE_DIR)

  # This runs the module with testing parameters:
  add_test(${PROJECT_NAME}_SampleTest
    ${RUNEXEC}
    --mi 50
    --cf -10
    --la 0.006
    --lsh 6
    --lcc test
    --numdir 26
    ${TEST_DATA_DIR}/dumbDWIData.nrrd
    ${TEST_DATA_DIR}/dumbSeeds.nii.gz
    ${OUTPUT_DIR}/finsler_distance_map.nii.gz
  )

  # This runs a program that checks the output and compares to the one expected:
  add_test(${PROJECT_NAME}_SampleCompare
    ${COMPEXEC}
    ${OUTPUT_DIR}/finsler_distance_map.nii.gz 
    ${TEST_DATA_DIR}/finsler_distance_map.nii.gz
  )

  add_test(${PROJECT_NAME}_SampleTest2
    ${RUNEXEC}
    --mi 50
    --cf -10
    --la 0.006
    --lsh 6
    --lcc test
    --numdir 26
    --useThreads
    ${TEST_DATA_DIR}/dumbDWIData.nrrd
    ${TEST_DATA_DIR}/dumbSeeds.nii.gz
    ${OUTPUT_DIR}/finsler_distance_map2.nii.gz
  )

  # This runs a program that checks the output and compares to the one expected:
  add_test(${PROJECT_NAME}_SampleCompare2
    ${COMPEXEC}
    ${OUTPUT_DIR}/finsler_distance_map2.nii.gz 
    ${TEST_DATA_DIR}/finsler_distance_map.nii.gz
  )

  add_test(${PROJECT_NAME}_SampleTest3
    ${RUNEXEC}
    --mi 50
    --cf -10
    --la 0.006
    --lsh 6
    --lcc test
    --numdir 26
    --useThreads
    --accel
    --accelIter 1
    ${TEST_DATA_DIR}/dumbDWIData.nrrd
    ${TEST_DATA_DIR}/dumbSeeds.nii.gz
    ${OUTPUT_DIR}/finsler_distance_map3.nii.gz
  )

  # This runs a program that checks the output and compares to the one expected:
  add_test(${PROJECT_NAME}_SampleCompare3
    ${COMPEXEC}
    ${OUTPUT_DIR}/finsler_distance_map3.nii.gz 
    ${TEST_DATA_DIR}/finsler_distance_map.nii.gz
  )


if( BUILD_BACKTRACING )
  add_test(${PROJECT_NAME}_SampleTest4
    ${RUNEXEC}
    --mi 50
    --cf -8
    --la 0.006
    --lsh 6
    --numdir 26
    --accel
    --accelIter 3
    --label 2
    --directionsImage ${OUTPUT_DIR}/sampleDirections.nrrd
    ${TEST_DATA_DIR}/sampleDWI.nrrd
    ${TEST_DATA_DIR}/sampleSeeds.nii.gz
    ${OUTPUT_DIR}/sampleFinslerCosts.nii.gz
  )
else( BUILD_BACKTRACING )
  add_test(${PROJECT_NAME}_SampleTest4
    ${RUNEXEC}
    --mi 50
    --cf -10
    --la 0.006
    --lsh 6
    --numdir 26
    --accel
    --accelIter 1
    --label 2
    ${TEST_DATA_DIR}/sampleDWI.nrrd
    ${TEST_DATA_DIR}/sampleSeeds.nii.gz
    ${OUTPUT_DIR}/sampleFinslerCosts.nii.gz
  )
endif( BUILD_BACKTRACING )

  # This runs a program that checks the output and compares to the one expected:
  add_test(${PROJECT_NAME}_SampleCompare4
    ${COMPEXEC}
    ${OUTPUT_DIR}/sampleFinslerCosts.nii.gz
    ${TEST_DATA_DIR}/sampleFinslerCosts.nii.gz
  )
  
  #This just removes the temporary files created
  add_test(${PROJECT_NAME}_CleanResults
    rm
    ${OUTPUT_DIR}/finsler_distance_map.nii.gz
    ${OUTPUT_DIR}/finsler_distance_map2.nii.gz
    ${OUTPUT_DIR}/finsler_distance_map3.nii.gz
    ${OUTPUT_DIR}/sampleFinslerCosts.nii.gz
  )

endif (BUILD_TESTING)
