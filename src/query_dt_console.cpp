//
// Copyright 2022 by Thomas Telkamp
//
// SPDX-License-Identifier: MIT
//

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <chrono>
#include <fstream>
#include <iostream>

#include <sys/time.h>

#include <assert.h>

// TCP/UDP
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>

int sockfd, connfd;

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

ssize_t readLine(int fd, void *buffer, size_t n) {

    ssize_t numRead;                    /* # of bytes fetched by last read() */
    size_t totRead;                     /* Total bytes read so far */
    char *buf;
    char ch;

    if (n <= 0 || buffer == NULL) {
        errno = EINVAL;
        return -1;
    }

    buf = (char *)buffer;                       /* No pointer arithmetic on "void *" */

    totRead = 0;
    for (;;) {
        numRead = read(fd, &ch, 1);

        if (numRead == -1) {
            if (errno == EINTR)         /* Interrupted --> restart read() */
                continue;
            else
                return -1;              /* Some other error */

        } else if (numRead == 0) {      /* EOF */
            if (totRead == 0)           /* No bytes read; return 0 */
                return 0;
            else                        /* Some bytes read; add '\0' */
                break;

        } else {                        /* 'numRead' must be 1 if we get here */
            if (totRead < n - 1) {      /* Discard > (n - 1) bytes */
                totRead++;
                *buf++ = ch;
            }

            if (ch == '\n')
                break;
        }
    }

    *buf = '\0';
    return totRead;
}

int main(int argc, char* argv[])
{
    // variables to be set by po
    std::string consolehost, format, pointing;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("console", po::value<std::string>(&consolehost)->default_value("console"), "console hostname")
        ("format", po::value<std::string>(&format)->default_value("text"), "output format (text, hms, filterbank)")
        ("pointing", po::value<std::string>(&pointing)->default_value("actual"), "RaDec pointing (actual, setpoint)")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("Get pointing data from the Dwingeloo Telescope console. %s") % desc << std::endl;
        std::cout << std::endl
                  << "This application gets RaDec and AzEl from the DT console.\n"
                  << std::endl;
        return ~0;
    }

    // bool bw_summary             = vm.count("progress") > 0;

    #define J2000_PORT 11030
    #define STAT_PORT  11042

    struct sockaddr_in servaddr, cli;

    // socket create and verification
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd == -1) {
        printf("socket creation failed...\n");
        exit(0);
    }

    int sockoptval = 1;
    setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &sockoptval, sizeof(int));
    sockoptval = 1;
    setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, &sockoptval, sizeof(int));

    /* don't wait when shutting down */
    linger lngr;
    lngr.l_onoff  = 1;
    lngr.l_linger = 0;
    setsockopt(sockfd, SOL_SOCKET, SO_LINGER, &lngr, sizeof(linger));

    bzero(&servaddr, sizeof(servaddr));

    struct hostent *host;

    if ((host = gethostbyname(consolehost.c_str())) == NULL)
    {
        printf("console hostname not known.\n");
        exit(0);
    }

    // assign IP, PORT
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = *(long *)(host->h_addr_list[0]);
    servaddr.sin_port = htons(J2000_PORT);

    // connect the client socket to server socket
    if (connect(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr))
        != 0) {
        printf("connection with the console failed.\n");
        exit(0);
    }
    else {
        // printf("connected to the console.\n");
    }

    char buffer[200];

    size_t len = readLine(sockfd, buffer, 200);


    buffer[len-1] = 0; // terminate string;

    std::string dt_radec(buffer);

    std::vector<std::string> radec_values;
    boost::split(radec_values, dt_radec, boost::is_any_of(" ,"));

    close(sockfd);

    // socket create and verification
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd == -1) {
        printf("socket creation failed...\n");
        exit(0);
    }

    setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &sockoptval, sizeof(int));
    setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, &sockoptval, sizeof(int));

    /* don't wait when shutting down */
    setsockopt(sockfd, SOL_SOCKET, SO_LINGER, &lngr, sizeof(linger));

    servaddr.sin_port = htons(STAT_PORT);

    // connect the client socket to server socket
    if (connect(sockfd, (struct sockaddr*)&servaddr, sizeof(servaddr))
        != 0) {
        printf("connection with the console (trace) failed.\n");
        exit(0);
    }
    // else {
    //     printf("connected to the console.\n");
    // }

    len = readLine(sockfd, buffer, 200);
    buffer[len-1] = 0; // terminate string;

    std::string dt_azel(buffer);

    std::vector<std::string> azel_values;
    boost::split(azel_values, dt_azel, boost::is_any_of(" ,"));

    char *ptr;
    double azimuth = 360.0*strtod(azel_values[0].c_str(), &ptr)/(2*M_PI) ;

    azimuth = fmod((azimuth+180.0),360.0);

    double elevation = 360.0*strtod(azel_values[1].c_str(), &ptr)/(2*M_PI) ;

    // OUTPUT

    uint8_t ra_index= 2;
    uint8_t dec_index = 3;

    if (pointing == "setpoint") {
        ra_index= 0;
        dec_index = 1;
    }

    if (format == "filterbank") {
        boost::erase_all(radec_values[ra_index], "d");
        boost::erase_all(radec_values[ra_index], "m");
        boost::erase_all(radec_values[ra_index], "h");
        boost::erase_all(radec_values[dec_index], "d");
        boost::erase_all(radec_values[dec_index], "m");
        boost::erase_all(radec_values[dec_index], "h");
        printf("%s,", radec_values[ra_index].c_str());
        printf("%s,", radec_values[dec_index].c_str());
        printf("%.3lf,", azimuth);
        printf("%.3lf\n", elevation);
    } else if (format == "hms") {
        printf("%s,", radec_values[ra_index].c_str());
        printf("%s,", radec_values[dec_index].c_str());
        printf("%.3lf,", azimuth);
        printf("%.3lf\n", elevation);
    } else if (format == "text") {
        printf("Right Ascension: %s\n", radec_values[ra_index].c_str());
        printf("Declination: %s\n", radec_values[dec_index].c_str());
        printf("Azimuth: %.2lf\n", azimuth);
        printf("Elevation: %.2lf\n", elevation);
    } else {
        printf ("unknown format\n");
    }

    /* clean up */
    // close the socket
    close(sockfd);

    // finished
    // std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
