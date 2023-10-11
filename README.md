# GIP

The Generic Integration Platform (GIP) is a proof of concept of my proposal to develop computational fluid dynamics(CFD) codes efficiently and highly portable.
This work is focused on solving the incompressible Navier-Stokes equations :
<p align="center"><img src="/imgs/NSequations.png" alt="drawing" width="300"/></p>


The main idea consists of reinterpreting the stencil-based operations as purely algebraic kernels. A special treatment is proposed for the non-linear operators to convert them into a concatenation of linear operators. The details of this implementation model are described in my article [Oyarzun2017](https://www.tandfonline.com/doi/abs/10.1080/10618562.2017.1390084).

The implementation structure is based on a set of abstract classes that provide a way of developing the computing node, parallel topology, and the core algebraic components of a CFD simulation.

- Parallel topology: describes how the data is distributed among the different computing nodes, and implements the communication scheme required to maintain consistency during the simulation.
- Nodes: specifies the computer architecture in which the code will run, this includes details about how many CPU cores, GPUs, and FPGAs can be used on each node.
- SpMV: the sparse matrix-vector multiplication is the dominant kernel in the CFD, in its specification special attention is needed in the sparsity pattern and storage format used.
- Vector operations: details the basic functions needed during the CFD simulation, mainly: axpy and dot.
- Timer: different architectures might implement the performance timers in different ways.

The linear solver and the CFD integrator are built as templates of the abstract classes, making them independent of the implementation details and improving their reusability on different computer architectures. These two classes are implemented once unless the simulation demands an other linear solver or specific integration scheme.

<p align="center"><img src="/imgs/structure.png" alt="drawing" width="300"/></p>

Whenever a new architecture wants to be tested, it is only necessary to implement the specification of the node, spmv, vector, and timer. This makes the codes highly modular and facilitates its portability across different platforms. 

The generation of the algebraic operators can be done in parallel or sequentially using any commercial or in-house CFD code.
In this repository, I have built an example of the implementation with specific implementations for pure MPI, MPI+OpenMP, MPI+CUDA, and MPI+OpenCL.

## Requirements
On your system you should have the following tools installed.

- MPI
- gcc
- cuda (mandatory if you plan to run with GPUs)

## Compiling in CPU mode

The compilation of the CPU version is as simple as running the script:
```
./compile_cpu.sh
```

After that you should get the executable file **./DrivenCavityCPU**

## Compiling in GPU mode
On GPU, you need to first build the GPU library
```
cd src/Specifications/Cuda/Libs/
./compileGpuLib.sh 
```
After that a file named **libgpu.a** should have been generated in the same folder.

Once you have it, you can run the compilation running the script:
```
./compile_gpu.sh
```
After that you should get the executable file **./DrivenCavityGPU**


## Execution 

The algebraic operators have been generated to distribute the computational domain in four parts. They are located in the InputFiles folder under the name C100K. If you plan to use your own operators just replace them.

So, the MPI must be executed with four processes like this:
```
mpirun -np 4 ./DrivenCavityCPU C100K
```
or, in the case running with GPUs:
```
mpirun -np 4 ./DrivenCavityGPU C100K
```

In both cases a correct execution should show this:
```
time 4.666345e-03 [dt::4.666345e-03] massr 1.728601e-06 ite 1
time 9.332691e-03 [dt::4.666345e-03] massr 8.792763e-07 ite 2
time 1.399904e-02 [dt::4.666345e-03] massr 9.428598e-07 ite 3
time 1.853262e-02 [dt::4.533584e-03] massr 8.464170e-07 ite 4
time 2.227019e-02 [dt::3.737568e-03] massr 6.869862e-07 ite 5
time 2.557171e-02 [dt::3.301522e-03] massr 5.971178e-07 ite 6
time 2.858014e-02 [dt::3.008431e-03] massr 5.321439e-07 ite 7
time 3.137420e-02 [dt::2.794062e-03] massr 4.830328e-07 ite 8
time 3.400284e-02 [dt::2.628633e-03] massr 4.448906e-07 ite 9
time 3.649912e-02 [dt::2.496283e-03] massr 4.138529e-07 ite 10
```


