#!/bin/bash

rmdir -rf build

mkdir build

cd build
cmake ..
cmake --build .
