# DIFI IQ tools
all: clients
clients: difi_to_file difi_fftmax difi_gnuplot difi_to_sigmf sigmf_to_difi difi_forwarder difi_spectrum difi_to_void
sdr: usrp_to_difi rfspace_to_difi rtlsdr_to_difi
gnuradio: difi_to_gnuradio

#INCLUDES = -I.
#LIBS = -L.

CFLAGS = -stdlib=libc++ -std=c++11
INCLUDES = -I. -I/opt/local/include -I../libvrt/include
LIBS = -L. -L../libvrt/build/ -L/opt/local/lib

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

difi_fftmax: difi_fftmax.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) -o difi_fftmax difi_fftmax.cpp \
		-lvrt -lzmq $(BOOSTLIBS) -lpthread -lfftw3

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

sigmf_to_difi: sigmf_to_difi.cpp
		g++ -O3 $(INCLUDES) $(LIBS) $(CFLAGS) sigmf_to_difi.cpp -o sigmf_to_difi \
		$(BOOSTLIBS) -lzmq -lvrt

convenience.o: convenience.c
		g++ -O3 -c $(INCLUDES) $(CFLAGS) -o convenience.o convenience.c

clean:
		$(RM) usrp_to_difi difi_to_file difi_fftmax difi_gnuplot difi_to_gnuradio difi_to_sigmf convenience.o rtlsdr_to_difi rfspace_to_difi difi_forwarder difi_to_void
