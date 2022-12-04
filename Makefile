# DIFI IQ tools
all: clients dt
clients: difi_to_file difi_fftmax difi_gnuplot difi_to_sigmf sigmf_to_difi difi_forwarder difi_spectrum difi_to_void control_difi difi_to_rtl_tcp difi_fftmax_quad difi_to_filterbank
sdr: usrp_to_difi rfspace_to_difi rtlsdr_to_difi
gnuradio: difi_to_gnuradio
gpu: difi_gpu_fftmax
dada: difi_to_dada
dt: query_dt_console

#INCLUDES = -I.
#LIBS = -L.

CFLAGS = -std=c++11
INCLUDES = -I. -I/opt/local/include -I../libvrt/include -I/opt/homebrew/include/
LIBS = -L. -L../libvrt/build/ -L/opt/local/lib -L/opt/homebrew/lib/

BOOSTLIBS = -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lboost_date_time
#BOOSTLIBS = -lboost_system-mt -lboost_program_options-mt -lboost_chrono-mt -lboost_filesystem-mt -lboost_thread-mt -lboost_date_time-mt

usrp_to_difi: usrp_to_difi.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) usrp_to_difi.cpp -o usrp_to_difi \
		-luhd -lpthread -lzmq -lvrt $(BOOSTLIBS)

difi_to_file: difi_to_file.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_file difi_to_file.cpp \
		-lvrt -lzmq -lpthread $(BOOSTLIBS)

difi_to_sigmf: difi_to_sigmf.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_sigmf difi_to_sigmf.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread

difi_to_gnuradio: difi_to_gnuradio.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_gnuradio difi_to_gnuradio.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lgnuradio-pmt

difi_to_void: difi_to_void.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_void difi_to_void.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

difi_to_difi_quad: difi_to_difi_quad.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_difi_quad difi_to_difi_quad.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

difi_to_rtl_tcp: difi_to_rtl_tcp.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_rtl_tcp difi_to_rtl_tcp.cpp \
		-lvrt -lzmq $(BOOSTLIBS)

difi_fftmax: difi_fftmax.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_fftmax difi_fftmax.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

difi_to_filterbank: difi_to_filterbank.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_filterbank difi_to_filterbank.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

difi_fftmax_quad: difi_fftmax_quad.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_fftmax_quad difi_fftmax_quad.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

# difi_pll_track: difi_pll_track.cpp
# 		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_pll_track difi_pll_track.cpp \
# 		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

difi_gnuplot: difi_gnuplot.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_gnuplot difi_gnuplot.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

difi_spectrum: difi_spectrum.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_spectrum difi_spectrum.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

difi_forwarder: difi_forwarder.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_forwarder difi_forwarder.cpp \
		-lzmq $(BOOSTLIBS)

rtlsdr_to_difi: convenience.o rtlsdr_to_difi.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) convenience.o rtlsdr_to_difi.cpp -o rtlsdr_to_difi \
		$(BOOSTLIBS) -lzmq -lvrt -lrtlsdr

rfspace_to_difi: rfspace_to_difi.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) rfspace_to_difi.cpp -o rfspace_to_difi \
		$(BOOSTLIBS) -lzmq -lvrt

query_dt_console: query_dt_console.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) query_dt_console.cpp -o query_dt_console \
		$(BOOSTLIBS)

sigmf_to_difi: sigmf_to_difi.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) sigmf_to_difi.cpp -o sigmf_to_difi \
		$(BOOSTLIBS) -lzmq -lvrt

control_difi: control_difi.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) control_difi.cpp -o control_difi \
		$(BOOSTLIBS) -lzmq -lvrt

difi_gpu_fftmax: difi_gpu_fftmax.cu
		nvcc -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_gpu_fftmax difi_gpu_fftmax.cu \
		$(BOOSTLIBS) -lzmq -lvrt -lcufft

difi_to_dada: difi_to_dada.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_to_dada difi_to_dada.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpsrdada \
		-I/home_local/camrasdemo/psrsoft/usr/include -L/home_local/camrasdemo/psrsoft/usr/lib

convenience.o: convenience.c
		g++ -O3 -c $(INCLUDES) $(CFLAGS) -o convenience.o convenience.c

clean:
		$(RM) usrp_to_difi difi_to_file difi_fftmax difi_gnuplot difi_to_gnuradio difi_to_sigmf convenience.o rtlsdr_to_difi rfspace_to_difi difi_forwarder difi_to_void difi_spectrum sigmf_to_difi difi_gpu_fftmax control_difi difi_to_dada difi_to_rtl_tcp difi_to_difi_quad difi_fftmax_quad difi_to_filterbank query_dt_console
