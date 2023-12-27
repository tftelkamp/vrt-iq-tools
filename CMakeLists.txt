cmake_minimum_required(VERSION 3.8)

project(VrtIqTools VERSION 1.0
                   DESCRIPTION "VRT IQ tools"
                   LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.50 REQUIRED COMPONENTS system program_options chrono filesystem thread date_time)
find_library(FFTW3_LIBRARY fftw3 REQUIRED)
find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h REQUIRED)

find_library(VRT_LIBRARY NAMES vrt REQUIRED)
find_path(VRT_INCLUDE_DIR NAMES "vrt/vrt_init.h" REQUIRED)

find_library(ZMQ_LIBRARY NAMES zmq REQUIRED)
find_path(ZMQ_INCLUDE_DIR NAMES "zmq.h" REQUIRED)

find_package(UHD)

if (NOT ${UHD_FOUND})
    message(WARNING "UHD not found, not building usrp_to_vrt")
endif()
find_library(GNURADIO_PMT_LIBRARY gnuradio-pmt QUIET)
if (NOT GNURADIO_PMT_LIBRARY)
    message(WARNING "Gnuradio-pmt library not found, not building vrt_to_gnuradio")
endif()

if (${UHD_FOUND})
    add_executable(usrp_to_vrt usrp_to_vrt.cpp)
    target_include_directories(usrp_to_vrt PRIVATE ${UHD_INCLUDE_DIRS})
    target_link_libraries(usrp_to_vrt PRIVATE ${UHD_LIBRARIES})
endif()

find_library(RTLSDR_LIBRARY rtlsdr)
find_path(RTLSDR_INCLUDE_DIR NAMES rtl-sdr.h)
if (NOT RTLSDR_LIBRARY)
    message(WARNING "RTL-SDR library library not found, not building rtlsdr_to_vrt")
else()
    set_source_files_properties(convenience.c PROPERTIES LANGUAGE CXX)
    add_executable(rtlsdr_to_vrt convenience.c rtlsdr_to_vrt.cpp)
    target_include_directories(rtlsdr_to_vrt PRIVATE ${RTLSDR_INCLUDE_DIR})
    target_link_libraries(rtlsdr_to_vrt PRIVATE ${RTLSDR_LIBRARY})
endif()

add_executable(vrt_to_sigmf vrt_to_sigmf.cpp)
add_executable(vrt_to_void vrt_to_void.cpp)
add_executable(vrt_to_udp vrt_to_udp.cpp)
add_executable(vrt_to_fifo vrt_to_fifo.cpp)
add_executable(vrt_pulsar vrt_pulsar.cpp)
add_executable(vrt_to_rtl_tcp vrt_to_rtl_tcp.cpp)
add_executable(vrt_fftmax vrt_fftmax.cpp)
add_executable(vrt_fftmax_quad vrt_fftmax_quad.cpp)
#add_executable(vrt_gpu_fftmax vrt_gpu_fftmax.cpp)
add_executable(vrt_gnuplot vrt_gnuplot.cpp)
add_executable(vrt_spectrum vrt_spectrum.cpp)
add_executable(vrt_forwarder vrt_forwarder.cpp)
add_executable(rfspace_to_vrt rfspace_to_vrt.cpp)
add_executable(query_dt_console query_dt_console.cpp)
add_executable(play_vrt play_vrt.cpp)
add_executable(control_vrt control_vrt.cpp)
add_executable(vrt_rffft vrt_rffft.cpp)
# TODO vrt_to_dada
# TODO rtlsdr_to_vrt

if (GNURADIO_PMT_LIBRARY)
    add_executable(vrt_to_gnuradio vrt_to_gnuradio.cpp)
    target_link_libraries(vrt_to_gnuradio PRIVATE ${GNURADIO_PMT_LIBRARY})
endif()

get_property(all_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)


foreach(target ${all_targets})
    target_include_directories(${target} PRIVATE ${FFTW3_INCLUDE_DIR})
    if (target STREQUAL "vrt_rffft")
        target_link_libraries(${target} PRIVATE ${FFTW3F_LIBRARY})
    else()
        target_link_libraries(${target} PRIVATE ${FFTW3_LIBRARY})
    endif()
    target_include_directories(${target} PRIVATE ${ZMQ_INCLUDE_DIR})
    target_link_libraries(${target} PRIVATE ${ZMQ_LIBRARY})
    target_include_directories(${target} PRIVATE ${VRT_INCLUDE_DIR})
    target_link_libraries(${target} PRIVATE ${VRT_LIBRARY})
    target_include_directories(${target} PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE ${Boost_LIBRARIES})
endforeach()


install()
