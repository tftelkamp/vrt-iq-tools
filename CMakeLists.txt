cmake_minimum_required(VERSION 3.8)

project(VrtIqTools VERSION 1.0
                   DESCRIPTION "VRT IQ tools"
                   LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.50 REQUIRED COMPONENTS system program_options chrono filesystem thread date_time)

find_library(VRT_LIBRARY NAMES vrt REQUIRED)
find_path(VRT_INCLUDE_DIR NAMES "vrt/vrt_init.h" REQUIRED)

find_library(ZMQ_LIBRARY NAMES zmq REQUIRED)
find_path(ZMQ_INCLUDE_DIR NAMES "zmq.h" REQUIRED)

find_package(UHD)

if (NOT ${UHD_FOUND})
    message(WARNING "UHD not found, not building usrp_to_vrt")
endif()
find_package(gnuradio-pmt QUIET)
if (NOT ${GNURADIO_FOUND})
    message(WARNING "GnuRadio not found, not building usrp_to_gnuradio")
endif()

if (${UHD_FOUND})
    add_executable(usrp_to_vrt usrp_to_vrt.cpp)
    target_include_directories(usrp_to_vrt PRIVATE ${UHD_INCLUDE_DIRS})
    target_link_libraries(usrp_to_vrt PRIVATE ${UHD_LIBRARIES})
endif()
add_executable(vrt_to_sigmf vrt_to_sigmf.cpp)
add_executable(vrt_to_gnuradio vrt_to_gnuradio.cpp)

get_property(all_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)


message(WARNING ${Boost_LIBRARIES})
foreach(target ${all_targets})
    target_include_directories(${target} PRIVATE ${ZMQ_INCLUDE_DIR})
    target_link_libraries(${target} PRIVATE ${ZMQ_LIBRARY})
    target_include_directories(${target} PRIVATE ${VRT_INCLUDE_DIR})
    target_link_libraries(${target} PRIVATE ${VRT_LIBRARY})
    target_include_directories(${target} PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE ${Boost_LIBRARIES})
endforeach()


install()
