# DIFI IQ tools
all: usrp_to_difi difi_to_file difi_fftmax difi_gnuplot difi_to_sigmf

INCLUDES = -I.
LIBS = -L.

usrp_to_difi: usrp_to_difi.cpp
		g++ -O3 $(INCLUDES) $(LIBS) usrp_to_difi.cpp -o usrp_to_difi \
		-luhd -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lpthread -lzmq -lvrt 

difi_to_file: difi_to_file.cpp
		g++ -O3 $(INCLUDES) $(LIBS) -o difi_to_file difi_to_file.cpp \
		-lvrt -lzmq -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lpthread

difi_to_sigmf: difi_to_sigmf.cpp
		g++ -O3 $(INCLUDES) $(LIBS) -o difi_to_sigmf difi_to_sigmf.cpp \
		-lvrt -lzmq -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lpthread

difi_to_gnuradio: difi_to_gnuradio.cpp
		g++ -O3 $(INCLUDES) $(LIBS) -o difi_to_gnuradio difi_to_gnuradio.cpp \
		-lvrt -lzmq -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lpthread -lgnuradio-pmt

difi_fftmax: difi_fftmax.cpp
		g++ -O3 $(INCLUDES) $(LIBS) -o difi_fftmax difi_fftmax.cpp \
		-lvrt -lzmq -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lpthread -lfftw3

difi_gnuplot: difi_gnuplot.cpp
		g++ -O3 $(INCLUDES) $(LIBS) -o difi_gnuplot difi_gnuplot.cpp \
		-lvrt -lzmq -lboost_system -lboost_program_options -lboost_chrono -lboost_filesystem -lboost_thread -lpthread -lfftw3
clean:
		$(RM) usrp_to_difi difi_to_file difi_fftmax difi_gnuplot difi_to_gnuradio
