#Cmake file for SuperCollider plugin
#some further Cmake hints at http://www.cs.swarthmore.edu/~adanner/tips/cmake.php
cmake_minimum_required (VERSION 2.6)
project (PolyPitch)

include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
#include_directories(${SC_PATH}/server/plugins) # for FFT_UGens.h; not needed from SC 3.5
include_directories(${SC_PATH}/external_libraries/libsndfile/)

set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

FIND_LIBRARY(Acc Accelerate)

include_directories(/usr/local/include/)	#samplerate.h
#include_directories(/opt/local/include/)	#sndfile.h

file(GLOB PolyPitchSources source/*cpp)
add_library(PolyPitch MODULE ${PolyPitchSources})

#if(APPLE)
#    set (PolyPitchLibs "-framework vecLib")
#else()
#    set (PolyPitchLibs ${FFTW3F_LIBRARY})
#endif()


#"/opt/local/lib/libsndfile.a" no need for libsndfile?
#64 bit version
target_link_libraries(PolyPitch "/usr/local/lib/libsamplerate.a" ${Acc} ${PolyPitchLibs})
#32 bit version
#target_link_libraries(PolyPitch "/Users/nickcollins/Desktop/tosort/libsamplerate32bit/libsamplerate.a" ${PolyPitchLibs})
