cmake_minimum_required(VERSION 3.14)

# ----- Basic properties -----
project (iMRS LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----- Linking with accelerate lib if on macos -----
IF(APPLE)
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
ENDIF(APPLE)

# ----- Find PZ. By default searches one directory up -----
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)
 

# ----- Check if PZ was installed with MKL -----
check_pz_opt(PZ_USING_MKL _PZ_HAS_MKL)
if(${_PZ_HAS_MKL} STREQUAL "OFF")
    message(FATAL_ERROR "\nERROR! NEOPZ WITHOUT MKL IS NOT SUPPORTED. PLEASE COMPILE NEOPZ WITH MKL")
endif()


file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})

# ----- To enable unit test using Catch2 UnitTest framework -----
option(BUILD_UNITTESTING "Whether the Catch2 UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
  include(cmake/add_unit_test.cmake)
  include(cmake/EnableCatch2.cmake)
  enable_catch2()
endif()

option(USING_DFN "Whether the DFN library will be linked in" ON)
if(USING_DFN)
    find_package(DFN REQUIRED HINTS "../../dfnResearch/dfnMesh_install/" ../dfnmesh_install/ ../../dfnmesh_install/ ${CMAKE_SOURCE_DIR}/dfnmesh_install/)
endif()



# ----- Adding subdirectory of iMRS library -----
add_subdirectory(Common)

# ----- Subdirectories that are most used -----
# add_subdirectory(BenchmarkCase1) Deprecated
add_subdirectory(OnlyFractures)
add_subdirectory(2Elements2Spaces2D)
# add_subdirectory(SimpleCaseTwoFractures)
add_subdirectory(SimpleCaseOneFracture)
add_subdirectory(CubicDomainNoFractures)
add_subdirectory(DFNIMRS)
add_subdirectory(Filling)
add_subdirectory(WANN)


# ----- Enables testing -----
if(BUILD_UNITTESTING)
    add_subdirectory(UnitTests)    
    ENABLE_TESTING()
    include(CTest)
endif()

# ----- Subdirectories that need fixing to new pz (with the crazy materials) -----
# add_subdirectory(SimpleCaseFracture)
# add_subdirectory(SimpleCaseTwoFractures)
# add_subdirectory(SimpleCaseThreeFractures)
# add_subdirectory(HDivSimple)

# ----- Subdirectories that need fixing and are seldom used -----
# add_subdirectory(MHMSimple)
# add_subdirectory(GeoMeshAndPropsGenerator)
# add_subdirectory(CGP)
# add_subdirectory(HdivCollapsed)
# add_subdirectory(Paper2D)