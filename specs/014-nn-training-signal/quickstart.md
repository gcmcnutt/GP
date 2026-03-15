# Quickstart: 014-nn-training-signal

## Prerequisites

```bash
sudo apt-get install -y libeigen3-dev libvtk9-dev g++ cmake gdb qtbase5-dev
# Note: Boost no longer required after refactoring
```

## Build

```bash
cd ~/GP/autoc
bash rebuild.sh
```

## Run NN Evolution (GA baseline)

```bash
cd ~/GP
./build/autoc
```

## Run NN Evolution (sep-CMA-ES)

Edit `autoc.ini`:
```ini
OptimizerType = sep-cma-es
PopulationSize = 50
NNMutationSigma = 0.3
NNSigmaFloor = 0.05
```

```bash
./build/autoc
```

## Run with Curriculum

Edit `autoc.ini`:
```ini
CurriculumEnabled = 1
CurriculumSchedule = 1:50,7:150,49:0
```

## Export to xiao-gp

```bash
./build/nnextractor -k nn-{timestamp} -o nn_weights.dat -i autoc.ini
./build/nn2cpp -i nn_weights.dat -o ~/xiao-gp/generated/nn_program_generated.cpp
cd ~/xiao-gp && pio run -e xiaoblesense_nn
```

## Run Tests

```bash
cd ~/GP/build
ctest --output-on-failure
```
