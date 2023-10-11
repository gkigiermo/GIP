# GIP

The Generic Integration Platform (GIP) is a proof of concept of my proposal to develop computational fluid dynamics(CFD) codes efficiently and highly portable.
This work is focused on solving the incompressible Navier-Stokes equations :
<p align="center"><img src="/imgs/NSequations.png" alt="drawing" width="300"/></p>


The main idea consists of reinterpreting the stencil-based operations as purely algebraic kernels. A special treatment is proposed for the non-linear operators to convert them into a concatenation of linear operators. The details of this implementation model are described in my article [Oyarzun2017](https://www.tandfonline.com/doi/abs/10.1080/10618562.2017.1390084).

The implementation structure is based on set of abstract classes that provides a way of developing the computing node, parallel topology, and the core algebraic componenents of a CFD simulation.

- Parallel topology: describes how the data is distributed among the different computing nodes, and implements the communication scheme required to maintain its consistency during the simulation.
- Nodes: specifies the computer architecture in which the code will run, this includes details about how many CPU-cores, GPUs, FPGAs can be used on each on each node.
- SpMV: the sparse matrix vector multiplication is the dominant kernel in the CFD, in its specification special attention is needed in the parsity pattern and storage format used.
- Vector operations: details the basic operations needed during the CFD simulation, mainly: axpy and dot.
- Timer: different architectures might implement the performance timers in different ways.

The linear solver and the CFD integrator are built as templates of the abstract classes, making them independent of the implementation details and improving its reusability on different computer architectures. This two classes are implemented once, unless the simulation demands a different linear solver or specific integration scheme.

<p align="center"><img src="/imgs/structure.png" alt="drawing" width="300"/></p>

Whenever a new architecture want to be tested, it is only necessary to implement the specification of the node, spmv, vector and timer. This makes the codes highly modular and facilitates its portability across different platforms. 

The generation of the algebraic operators can be done in parallel or sequentially using any comercial or in-house CFD code.
In this repository, I have built an example of the implementation with specific implementations for : pure MPI, MPI+OpenMP, MPI+CUDA and MPI+OpenCL.
