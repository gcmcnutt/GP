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

Also need to set Xvfb for testing and headless operations
```
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
./your_opengl_application
```
