#!/bin/sh
g++ --std=c++11 inc/model.cpp inc/k-means.cpp arkr.cpp -I ../boost_1_57_0/ -O3 -o arkr
