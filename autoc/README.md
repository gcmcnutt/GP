# Building autoc

# on u22
sudo apt-get update
sudo apt-get install -y libboost-all-dev libeigen3-dev libvtk9-dev xvfb g++ cmake gdb qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools

# install asw cli
# establish credentials (aws configure)
# install aws-sdk-vpp https://github.com/aws/aws-sdk-cpp


...

- First make sure top level GP is built
- Generate local makefile (vscode can auto do this)
```
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
