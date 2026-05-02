cmake_minimum_required(VERSION 3.17)

include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

project(
  VrtIqTools
  VERSION 1.0
  DESCRIPTION "VRT IQ tools"
  LANGUAGES CXX)

set(MIMICPP_VERSION v9.3.0)
set(CATCH2_VERSION v3.8.1)

# 20 or higher is needed for mimicpp
set(CMAKE_CXX_STANDARD 20)

# To use clangd as LSP (for Neovim / Vim etc)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compile time options
OPTION(VRT_IQ_TOOLS_TESTING OFF)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(Boost 1.50 REQUIRED COMPONENTS program_options chrono filesystem
                                            thread date_time)
find_library(FFTW3_LIBRARY fftw3 REQUIRED)
find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
find_library(FFTW3_THREADS_LIBRARY fftw3_threads REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h REQUIRED)

# Fetch the submodule if not found
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    message(STATUS "Fetching libvrt submodule...")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT
    )
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "Failed to fetch libvrt submodule. Git returned error code ${GIT_SUBMOD_RESULT}")
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT GIT_BRANCH)
      set(GIT_BRANCH "unknown")
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT GIT_COMMIT)
      set(GIT_COMMIT "unknown")
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --format=%ci
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DATE
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT GIT_DATE)
      set(GIT_DATE "unknown")
    endif()

    add_compile_definitions(GIT_BRANCH="${GIT_BRANCH}" GIT_COMMIT="${GIT_COMMIT}" GIT_DATE="${GIT_DATE}")
else()
    message(FATAL_ERROR "Git is required to fetch libvrt submodule. Please make sure Git is installed.")
endif()
add_subdirectory(libvrt)

include_directories(${CMAKE_SOURCE_DIR}/include)

find_library(ZMQ_LIBRARY NAMES zmq REQUIRED)
find_path(ZMQ_INCLUDE_DIR NAMES "zmq.h" REQUIRED)

find_package(UHD QUIET)
if(${UHD_FOUND})
  add_executable(usrp_to_vrt src/usrp_to_vrt.cpp)
  target_include_directories(usrp_to_vrt PRIVATE ${UHD_INCLUDE_DIRS})
  target_link_libraries(usrp_to_vrt PRIVATE ${UHD_LIBRARIES})
else()
  message(WARNING "UHD not found, not building usrp_to_vrt")
endif()

find_package(HackRF)
if(${HackRF_FOUND})
  add_executable(hackrf_to_vrt hackrf_to_vrt.cpp)
  target_link_libraries(hackrf_to_vrt PRIVATE HackRF::hackrf_static)
else()
  message(WARNING "HackRF not found, not building hackrf_to_vrt")
endif()

find_library(RTLSDR_LIBRARY rtlsdr)
if(RTLSDR_LIBRARY)
  find_path(RTLSDR_INCLUDE_DIR NAMES rtl-sdr.h REQUIRED)
  set_source_files_properties(src/convenience.c PROPERTIES LANGUAGE CXX)
  add_executable(rtlsdr_to_vrt src/convenience.c src/rtlsdr_to_vrt.cpp)
  target_include_directories(rtlsdr_to_vrt PRIVATE ${RTLSDR_INCLUDE_DIR})
  target_link_libraries(rtlsdr_to_vrt PRIVATE ${RTLSDR_LIBRARY})
else()
  message(
    WARNING "RTL-SDR library library not found, not building rtlsdr_to_vrt")
endif()

find_library(AIRSPY_LIBRARY airspy)
if(AIRSPY_LIBRARY)
  find_path(AIRSPY_INCLUDE_DIR NAMES libairspy/airspy.h REQUIRED)
  add_executable(airspy_to_vrt src/airspy_to_vrt.cpp)
  target_include_directories(airspy_to_vrt PRIVATE ${AIRSPY_INCLUDE_DIR})
  target_link_libraries(airspy_to_vrt PRIVATE ${AIRSPY_LIBRARY})
else()
  message(
    WARNING "Airspy library not found, not building airspy_to_vrt")
endif()

find_library(IIO_LIBRARY iio)
find_library(AD9361_LIBRARY ad9361)
if(IIO_LIBRARY AND AD9361_LIBRARY)
  find_path(IIO_INCLUDE_DIR NAMES iio.h REQUIRED)
  find_path(AD9361_INCLUDE_DIR NAMES ad9361.h REQUIRED)
  add_executable(iio_to_vrt src/iio_to_vrt.cpp)
  target_include_directories(iio_to_vrt PRIVATE ${IIO_INCLUDE_DIR} ${AD9361_INCLUDE_DIR})
  target_link_libraries(iio_to_vrt PRIVATE ${IIO_LIBRARY} ${AD9361_LIBRARY})
else()
  message(
    WARNING "libiio and libad9361 not both found, building without iio_to_vrt (pluto)")
endif()

find_package(CUDAToolkit)
if(NOT CUDAToolkit_FOUND)
  message(WARNING "CUDA not found, building without GPU support")
else()
  enable_language(CUDA)
  find_library(
    CUFFT_LIBRARY
    NAMES cufft
    HINTS ${CUDAToolkit_LIBRARY_DIR})
  add_executable(vrt_gpu_fftmax src/vrt_gpu_fftmax.cu)
  set_property(TARGET vrt_gpu_fftmax PROPERTY CUDA_ARCHITECTURES OFF)
  target_link_libraries(vrt_gpu_fftmax PRIVATE ${CUFFT_LIBRARY})
  add_executable(vrt_gpu_channelizer src/vrt_gpu_channelizer.cu)
  set_property(TARGET vrt_gpu_channelizer PROPERTY CUDA_ARCHITECTURES OFF)
  target_link_libraries(vrt_gpu_channelizer PRIVATE ${CUFFT_LIBRARY})
endif()

add_executable(vrt_to_sigmf src/vrt_to_sigmf.cpp)
add_executable(vrt_version src/vrt_version.cpp)
add_executable(sigmf_to_vrt src/sigmf_to_vrt.cpp)
add_executable(vrt_to_void src/vrt_to_void.cpp)
add_executable(vrt_to_stdout src/vrt_to_stdout.cpp)
add_executable(vrt_to_udp src/vrt_to_udp.cpp)
add_executable(vrt_to_fifo src/vrt_to_fifo.cpp)
add_executable(vrt_to_filterbank src/vrt_to_filterbank.cpp)
add_executable(vrt_pulsar src/vrt_pulsar.cpp)
add_executable(vrt_to_rtl_tcp src/vrt_to_rtl_tcp.cpp)
add_executable(vrt_fftmax src/vrt_fftmax.cpp)
add_executable(vrt_fftmax_quad src/vrt_fftmax_quad.cpp)
add_executable(vrt_spectrum src/vrt_spectrum.cpp)
add_executable(vrt_forwarder src/vrt_forwarder.cpp)
add_executable(rfspace_to_vrt src/rfspace_to_vrt.cpp)
add_executable(query_dt_console src/query_dt_console.cpp)
add_executable(play_vrt src/play_vrt.cpp)
add_executable(control_vrt src/control_vrt.cpp)
add_executable(vrt_rffft src/vrt_rffft.cpp)
add_executable(vrt_metadata src/vrt_metadata.cpp)
add_executable(vrt_tuner src/vrt_tuner.cpp)
add_executable(vrt_merge src/vrt_merge.cpp)
add_executable(vrt_buffer src/vrt_buffer.cpp)
add_executable(vrt_correlate src/vrt_correlate.cpp)

find_library(GNURADIO_PMT_LIBRARY gnuradio-pmt QUIET)
if(GNURADIO_PMT_LIBRARY)
  add_executable(vrt_to_gnuradio src/vrt_to_gnuradio.cpp)
  target_link_libraries(vrt_to_gnuradio PRIVATE ${GNURADIO_PMT_LIBRARY})
else()
  message(
    WARNING "Gnuradio-pmt library not found, not building vrt_to_gnuradio")
endif()

find_library(PSRDADA_LIBRARY psrdada QUIET)
if(PSRDADA_LIBRARY)
  find_path(PSRDADA_INCLUDE_DIR NAMES dada_def.h REQUIRED)
  add_executable(vrt_to_dada src/vrt_to_dada.cpp)
  target_include_directories(vrt_to_dada PRIVATE ${PSRDADA_INCLUDE_DIR})
  target_link_libraries(vrt_to_dada PRIVATE ${PSRDADA_LIBRARY})
else()
  message(WARNING "PSRDADA not found, not building vrt_to_dada")
endif()

get_property(
  all_targets
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY BUILDSYSTEM_TARGETS)

foreach(target ${all_targets})
  target_include_directories(${target} PRIVATE ${FFTW3_INCLUDE_DIR})
  if(target STREQUAL "vrt_rffft")
    target_link_libraries(${target} PRIVATE ${FFTW3F_LIBRARY})
  else()
    target_link_libraries(${target} PRIVATE ${FFTW3_LIBRARY})
  endif()
  if(target STREQUAL "vrt_to_filterbank")
    target_link_libraries(${target} PRIVATE ${FFTW3_THREADS_LIBRARY})
  endif()
  target_include_directories(${target} PRIVATE ${ZMQ_INCLUDE_DIR})
  target_link_libraries(${target} PRIVATE ${ZMQ_LIBRARY})
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/libvrt/include)
  target_link_libraries(${target} PRIVATE vrt)
  target_include_directories(${target} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(${target} PRIVATE ${Boost_LIBRARIES})
endforeach()

install(TARGETS ${all_targets})


if(VRT_IQ_TOOLS_TESTING)
    FetchContent_Declare(mimicpp
        GIT_REPOSITORY https://github.com/DNKpp/mimicpp
        VERSION ${MIMICPP_VERSION}
        GIT_TAG ${MIMICPP_VERSION}
    )
    FetchContent_MakeAvailable(mimicpp)

    find_package(Catch2 3 QUIET)
    if (NOT Catch2_FOUND)
        message("Catch2 not found, downloading version: ${CATCH2_VERSION} ...")
        FetchContent_Declare(Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG ${CATCH2_VERSION}
        )
        FetchContent_MakeAvailable(Catch2)
    endif ()
endif()