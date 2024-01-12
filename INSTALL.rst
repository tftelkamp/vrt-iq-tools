Intall vrt-iq-tools
===================

Dependencies
************

C++ Boost libraries
-------------------

Install the boost development packages using the package manager of the 
destribution:

- boost_system 
- boost_program_options 
- boost_chrono
- boost_filesystem
- boost_thread
- boost_date_time

Libraries
---------

The libvrt is available on github and 
can be cloned using::

    git clone https://github.com/tftelkamp/libvrt

When using the rtlsdr dongle, install:

- rtl-sdr-devel

The Ettus SDR requires:

- uhd-devel


Build vrt-iq-tools
******************

Use make and tab completion to list all the make targets. To install the vrt-iq-tools, use::

    make
    
For SDR targets run::

    make sdr
    
And if you use thr RTLSDR device, you can use::

    make rtlsdr

The tools are installed in /usr/local/bin using::

    sudo make install

 

