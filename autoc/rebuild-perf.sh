#! /bin/bash
set -ex

# no parallelism, no special knobs, force hard clean 

# rebuild from clean
cd ~/GP/autoc
rm -rf build
mkdir build
cd build
cmake -DPERFORMANCE_BUILD=ON ..
make

# same for crrcsim -- perf mode
cd ~/crsim/crrcsim-0.9.13
rm -rf build
mkdir build
cd build
cmake -DPERFORMANCE_BUILD=ON ..
make
