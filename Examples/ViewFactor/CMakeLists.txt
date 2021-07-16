cmake_minimum_required(VERSION 3.4)

project("ViewFactor")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

## Enable Clang sanitizer for debug builds
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" CACHE STRING "")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" CACHE STRING "")
endif()

# set default build type
SET(DEFAULT_BUILD_TYPE "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# automatically append standard path for ViennaHRLEConfig.cmake
# might need to be changed for different install directories
# not necessary when installed to /usr/local (not recommended)
SET(ViennaToolsDir "/home/filipov/Software/ViennaTools/")
SET(ViennaHRLE_DIR "${ViennaToolsDir}/ViennaHRLE/install/ViennaHRLE/lib/cmake/ViennaHRLE")
SET(ViennaLS_DIR "${ViennaToolsDir}/ViennaLS/install/ViennaLS/lib/cmake/ViennaLS")

find_package(ViennaHRLE REQUIRED)
find_package(ViennaLS REQUIRED)

#option(BUILD_TEST "Build unit test." OFF)
#if(BUILD_TEST)
#  add_executable(Test Test/Test.cpp)
#  target_include_directories(Test PUBLIC ${VIENNALS_INCLUDE_DIRS})
#  target_link_libraries(Test PRIVATE ${VIENNALS_LIBRARIES})
#endif()

set(ViewFactor "ViewFactor")
add_executable(${ViewFactor} ${ViewFactor}.cpp)
target_include_directories(${ViewFactor} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${ViewFactor} PRIVATE ${VIENNALS_LIBRARIES})

