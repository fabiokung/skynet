#!/bin/bash

SKYNET_INSTALL=/path/to/SkyNet/install

g++ -std=c++11 calc_error.cpp -I/usr/include/hdf5/serial -I${SKYNET_INSTALL}/include/ -L${SKYNET_INSTALL}/lib/ -lhdf5_serial -lhdf5_cpp -lboost_filesystem -lboost_system -lSkyNet_static -o calc_error
