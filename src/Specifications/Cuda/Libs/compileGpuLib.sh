#!/bin/bash
echo "Compiling ..."
nvcc -arch sm_70 -c -o libgpu.a libgpu.cu --compiler-options '-fPIC'
echo "Installing library"
cp GIP_gpu.h ../Headers/.
cp libgpu.a ../Libs/.
echo "Ready for Use!!"
