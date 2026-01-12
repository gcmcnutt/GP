#! /bin/bash
set -ex

# no parallelism, no special knobs, force hard clean 

# rebuild from clean
cd ~/GP
make clean
rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../autoc
cd ..
make

# same for crrcsim
cd ~/crsim/crrcsim-0.9.13
rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
