# Building autoc

need libboost-all-dev
need Eigen
need vtk
...

- First make sure top level GP is built
- Generate local makefile
```
cd autoc
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```
