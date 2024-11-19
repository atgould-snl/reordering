#!/bin/bash
mkdir build
cd build
TRILBUILD=/scratch/malphil/trilinos-devel/trilinos-builds
TRILINOS=/scratch/malphil/trilinos-devel/Trilinos
TRILINSTALL=$TRILBUILD/installs/ifpack2
cmake \
  -DTrilinos_DIR=$TRILINSTALL/lib64/cmake/Trilinos \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -DCMAKE_CXX_EXTENSIONS=OFF \
  ..

make -j 12
