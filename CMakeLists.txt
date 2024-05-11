cmake_minimum_required(VERSION 3.17)

project(
  VrtIqTools
  VERSION 1.0
  DESCRIPTION "VRT IQ tools"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.50 REQUIRED COMPONENTS system program_options chrono
                                            filesystem thread date_time)
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
else()
    message(FATAL_ERROR "Git is required to fetch libvrt submodule. Please make sure Git is installed.")
endif()
add_subdirectory(libvrt)

find_library(ZMQ_LIBRARY NAMES zmq REQUIRED)
find_path(ZMQ_INCLUDE_DIR NAMES "zmq.h" REQUIRED)

find_package(UHD QUIET)
if(${UHD_FOUND})
  add_executable(usrp_to_vrt usrp_to_vrt.cpp)
  target_include_directories(usrp_to_vrt PRIVATE ${UHD_INCLUDE_DIRS})
  target_link_libraries(usrp_to_vrt PRIVATE ${UHD_LIBRARIES})
else()
  message(WARNING "UHD not found, not building usrp_to_vrt")
endif()

find_library(RTLSDR_LIBRARY rtlsdr)
if(RTLSDR_LIBRARY)
  find_path(RTLSDR_INCLUDE_DIR NAMES rtl-sdr.h REQUIRED)
  set_source_files_properties(convenience.c PROPERTIES LANGUAGE CXX)
  add_executable(rtlsdr_to_vrt convenience.c rtlsdr_to_vrt.cpp)
  target_include_directories(rtlsdr_to_vrt PRIVATE ${RTLSDR_INCLUDE_DIR})
  target_link_libraries(rtlsdr_to_vrt PRIVATE ${RTLSDR_LIBRARY})
else()
  message(
    WARNING "RTL-SDR library library not found, not building rtlsdr_to_vrt")
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
  add_executable(vrt_gpu_fftmax vrt_gpu_fftmax.cu)
  set_property(TARGET vrt_gpu_fftmax PROPERTY CUDA_ARCHITECTURES OFF)
  target_link_libraries(vrt_gpu_fftmax PRIVATE ${CUFFT_LIBRARY})
endif()

add_executable(vrt_to_sigmf vrt_to_sigmf.cpp)
add_executable(sigmf_to_vrt sigmf_to_vrt.cpp)
add_executable(vrt_to_void vrt_to_void.cpp)
add_executable(vrt_to_stdout vrt_to_stdout.cpp)
add_executable(vrt_to_udp vrt_to_udp.cpp)
add_executable(vrt_to_fifo vrt_to_fifo.cpp)
add_executable(vrt_to_filterbank vrt_to_filterbank.cpp)
add_executable(vrt_pulsar vrt_pulsar.cpp)
add_executable(vrt_to_rtl_tcp vrt_to_rtl_tcp.cpp)
add_executable(vrt_fftmax vrt_fftmax.cpp)
add_executable(vrt_fftmax_quad vrt_fftmax_quad.cpp)
add_executable(vrt_gnuplot vrt_gnuplot.cpp)
add_executable(vrt_spectrum vrt_spectrum.cpp)
add_executable(vrt_forwarder vrt_forwarder.cpp)
add_executable(rfspace_to_vrt rfspace_to_vrt.cpp)
add_executable(query_dt_console query_dt_console.cpp)
add_executable(play_vrt play_vrt.cpp)
add_executable(control_vrt control_vrt.cpp)
add_executable(vrt_rffft vrt_rffft.cpp)
add_executable(vrt_metadata vrt_metadata.cpp)

find_library(GNURADIO_PMT_LIBRARY gnuradio-pmt QUIET)
if(GNURADIO_PMT_LIBRARY)
  add_executable(vrt_to_gnuradio vrt_to_gnuradio.cpp)
  target_link_libraries(vrt_to_gnuradio PRIVATE ${GNURADIO_PMT_LIBRARY})
else()
  message(
    WARNING "Gnuradio-pmt library not found, not building vrt_to_gnuradio")
endif()

find_library(PSRDADA_LIBRARY psrdada QUIET)
if(PSRDADA_LIBRARY)
  find_path(PSRDADA_INCLUDE_DIR NAMES dada_def.h REQUIRED)
  add_executable(vrt_to_dada vrt_to_dada.cpp)
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
