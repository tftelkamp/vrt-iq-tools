# PREFIX is environment variable, but if it is not set, then set default value
ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

CXX=g++

# VRT IQ tools
all: clients dt
clients: vrt_fftmax vrt_gnuplot vrt_to_sigmf sigmf_to_vrt play_vrt vrt_forwarder vrt_spectrum vrt_to_void control_vrt vrt_to_rtl_tcp vrt_fftmax_quad vrt_to_filterbank vrt_to_fifo vrt_pulsar vrt_to_udp vrt_metadata vrt_to_stdout
sdr: usrp_to_vrt rfspace_to_vrt rtlsdr_to_vrt
gnuradio: vrt_to_gnuradio
gpu: vrt_gpu_fftmax
dada: vrt_to_dada
dt: query_dt_console
strf: vrt_rffft

#INCLUDES = -I.
#LIBS = -L.

CFLAGS = -std=c++11
INCLUDES = -I. -I/opt/local/include -I../libvrt/include -I/opt/homebrew/include/
LIBS = -L. -L../libvrt/build/ -L/usr/local/lib -L/opt/local/lib -L/opt/homebrew/lib/

BOOSTLIBS = -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lboost_date_time
# BOOSTLIBS = -lboost_system-mt -lboost_program_options-mt -lboost_chrono-mt -lboost_filesystem-mt -lboost_thread-mt -lboost_date_time-mt

usrp_to_vrt: usrp_to_vrt.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) usrp_to_vrt.cpp -o usrp_to_vrt \
		-luhd -lpthread -lzmq -lvrt $(BOOSTLIBS)

vrt_to_sigmf: vrt_to_sigmf.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_sigmf vrt_to_sigmf.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread

vrt_to_gnuradio: vrt_to_gnuradio.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_gnuradio vrt_to_gnuradio.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lgnuradio-pmt

vrt_to_void: vrt_to_void.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_void vrt_to_void.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

vrt_to_stdout: vrt_to_stdout.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_stdout vrt_to_stdout.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

vrt_to_udp: vrt_to_udp.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_udp vrt_to_udp.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

vrt_to_fifo: vrt_to_fifo.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_fifo vrt_to_fifo.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

vrt_to_vrt_quad: vrt_to_vrt_quad.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_vrt_quad vrt_to_vrt_quad.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

vrt_to_rtl_tcp: vrt_to_rtl_tcp.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_rtl_tcp vrt_to_rtl_tcp.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

vrt_fftmax: vrt_fftmax.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_fftmax vrt_fftmax.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

vrt_pulsar: vrt_pulsar.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_pulsar vrt_pulsar.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

vrt_to_filterbank: vrt_to_filterbank.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_filterbank vrt_to_filterbank.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3 -lfftw3_threads

vrt_fftmax_quad: vrt_fftmax_quad.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_fftmax_quad vrt_fftmax_quad.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

vrt_gnuplot: vrt_gnuplot.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_gnuplot vrt_gnuplot.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

vrt_spectrum: vrt_spectrum.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_spectrum vrt_spectrum.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

vrt_metadata: vrt_metadata.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_metadata vrt_metadata.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

vrt_forwarder: vrt_forwarder.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_forwarder vrt_forwarder.cpp \
		-lzmq $(BOOSTLIBS)

rtlsdr_to_vrt: convenience.o rtlsdr_to_vrt.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) convenience.o rtlsdr_to_vrt.cpp -o rtlsdr_to_vrt \
		$(BOOSTLIBS) -lzmq -lvrt -lrtlsdr

rfspace_to_vrt: rfspace_to_vrt.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) rfspace_to_vrt.cpp -o rfspace_to_vrt \
		$(BOOSTLIBS) -lzmq -lvrt

query_dt_console: query_dt_console.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) query_dt_console.cpp -o query_dt_console \
		$(BOOSTLIBS)

sigmf_to_vrt: sigmf_to_vrt.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) sigmf_to_vrt.cpp -o sigmf_to_vrt \
		$(BOOSTLIBS) -lzmq -lvrt

play_vrt: play_vrt.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) play_vrt.cpp -o play_vrt \
		$(BOOSTLIBS) -lzmq -lvrt

control_vrt: control_vrt.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) control_vrt.cpp -o control_vrt \
		$(BOOSTLIBS) -lzmq -lvrt

vrt_gpu_fftmax: vrt_gpu_fftmax.cu
		nvcc -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_gpu_fftmax vrt_gpu_fftmax.cu \
		$(BOOSTLIBS) -lzmq -lvrt -lcufft

vrt_to_dada: vrt_to_dada.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o vrt_to_dada vrt_to_dada.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpsrdada \
		-I/home_local/camrasdemo/psrsoft/usr/include -L/home_local/camrasdemo/psrsoft/usr/lib

vrt_rffft: vrt_rffft.cpp
		${CXX} -O3 $(INCLUDES) $(LIBS) $(CFLAGS) vrt_rffft.cpp -o vrt_rffft \
		$(BOOSTLIBS) -lzmq -lvrt -lfftw3f

convenience.o: convenience.c
		${CXX} -O3 -c $(INCLUDES) $(CFLAGS) -o convenience.o convenience.c

install: all
		install -m 755 vrt_fftmax        $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_gnuplot       $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_to_sigmf      $(DESTDIR)$(PREFIX)/bin/
		install -m 755 sigmf_to_vrt      $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_forwarder     $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_spectrum      $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_to_void       $(DESTDIR)$(PREFIX)/bin/
		install -m 755 control_vrt       $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_to_rtl_tcp    $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_fftmax_quad   $(DESTDIR)$(PREFIX)/bin/
		install -m 755 vrt_to_filterbank $(DESTDIR)$(PREFIX)/bin/
		install -m 755 query_dt_console   $(DESTDIR)$(PREFIX)/bin/

clean:
		$(RM) usrp_to_vrt vrt_fftmax vrt_gnuplot vrt_to_gnuradio vrt_to_sigmf convenience.o rtlsdr_to_vrt rfspace_to_vrt vrt_forwarder vrt_to_void vrt_spectrum sigmf_to_vrt play_vrt vrt_gpu_fftmax control_vrt vrt_to_dada vrt_to_rtl_tcp vrt_to_vrt_quad vrt_fftmax_quad vrt_to_filterbank query_dt_console vrt_rffft vrt_to_fifo vrt_pulsar vrt_to_udp vrt_metadata vrt_to_stdout
