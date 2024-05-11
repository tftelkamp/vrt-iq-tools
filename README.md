# vrt-iq-tools

vrt-iq-tools is a collection of programs used for streaming IQ data using (a slightly modified variant of) [libvrt](https://github.com/ember91/libvrt), which in turn implements part of the [ANSI/VITA 49.0 Radio Transport](https://www.vita.com/resources/Documents/Articles/IEEE%20CAS%20Mag%202012.pdf) (VRT) standard. Most of the tools use VRT streams over ZMQ sockets.

## Installation

### Dependencies

This package needs much of libboost, libzmq, libfftw3 and libfftw3f. On a clean Ubuntu 22.04 system, these can be installed with:


```bash
apt install libboost-chrono-dev libboost-date-time-dev libboost-filesystem-dev libboost-program-options-dev libboost-system-dev libboost-thread-dev libzmq3-dev libfftw3-dev libfftw3-single3
```

Optional dependencies:

* `librtlsdr-dev` for RTL-SDR support
* `libuhd-dev` for UHD (USRP) support
* `libgnuradio-pmt3.10.1` for gnuradio support
* `nvidia-cuda-toolkit-gcc` for GPU acceleration of some tools
* `libpsrdada` for DADA support, see [installation instructions](https://psrdada.sourceforge.net/download.shtml)

### Compilation

```bash
git clone https://github.com/tftelkamp/vrt-iq-tools
cd vrt-iq-tools
mkdir build
cd build
cmake ..
make install -j
```

## Usage

### Creating a VRT stream from an SDR or file:

* `usrp_to_vrt`: Create VRT stream from an [Ettus USRP](https://www.ettus.com/products/) SDR device (e.g. B210).
* `rfspace_to_vrt`: Create VRT stream from [RFSpace](https://http://www.rfspace.com) SDR device.
* `rtlsdr_to_vrt`: Create VRT stream from [RTL-SDR](https://www.rtl-sdr.com/) device.
* `sigmf_to_vrt`: Create VRT stream from [SigMF](https://sigmf.org) recording, or with `--vrt` from a VRT recording.
* `play_vrt`: Create VRT stream from [SigMF](https://sigmf.org) recording, intended for transmitting.

### Clients:

* `vrt_to_sigmf`: Store IQ and metadata as [SigMF](https://sigmf.org) recording, or with `--vrt` as raw VRT.
* `vrt_spectrum`: Create spectra, store in CSV or ECSV format (compatible with [Astropy](https://astropy.org)). With `--gnuplot`, output can be piped to Gnuplot.
* `vrt_to_filterbank`: Create spectra, store in [sigproc](https://sigproc.sourceforge.net/) filterbank format.
* `vrt_rffft`: Create spectra and store in [STRF](https://github.com/cbassa/strf) format.
* `vrt_fftmax`: Create spectra, store only the frequency of the bin with the maximum. Used for Doppler tracking.
* `vrt_fftmax_quad`: Same as `vrt_fftmax` but used for modulated signals.
* `vrt_pulsar`: Channelize, dedisperse and fold pulsar data.

### Converting to other stream types:
* `vrt_to_stdout`: Stream IQ to standard output. Useful for streaming to [PhantomSDR](https://github.com/PhantomSDR/PhantomSDR).
* `vrt_to_rtl_tcp`: Stream as 8-bit RTL-TCP stream.
* `vrt_to_gnuradio`: Stream IQ to ZeroMQ socket to be used in [GNURadio](https://www.gnuradio.org).
* `vrt_to_fifo`: Write IQ to a fifo buffer.
* `vrt_to_udp`: Stream IQ as fc32 UDP packets.
* `vrt_to_dada`: Stream to a [psrdada](https://psrdada.sourceforge.net/) ring buffer, to be used with pulsar software such as [dspsr](https://dspsr.sourceforge.net/).

### Miscellaneous:
* `vrt_metadata`: Print metadata of a VRT stream.
* `vrt_forwarder`: Forward ZMQ stream.
* `vrt_to_void`: Template for new clients.
* `control_vrt`: Control devices, e.g. to set gain or frequency.

## License

[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)