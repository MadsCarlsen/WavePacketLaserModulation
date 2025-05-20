cmake_minimum_required(VERSION 3.24)
project(laser_modulation)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -march=native")  # -fopt-info-vec-optimized -fopt-info-vec-missed 

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find Boost 
find_package(Boost REQUIRED)

# Find FFTW3 (double-precision)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

add_executable(laser_modulation main.cpp
        EvanescentField.cpp
        EvanescentField.h
        HelperFunctions.cpp
        HelperFunctions.h
        WavePacket.cpp
        WavePacket.h
        FourierTransformation.cpp
        FourierTransformation.h
        Grid2D.h
        Spline2D.h
        WavePacketModulation.cpp
        WavePacketModulation.h
        BilinearInterpolatorRegular.h
        BicubicInterpolator.h
)

message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")

# Include directories 
target_include_directories(laser_modulation PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
)

target_link_libraries(laser_modulation PRIVATE
        PkgConfig::FFTW
        OpenMP::OpenMP_CXX
)
