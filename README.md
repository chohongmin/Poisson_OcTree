# Poisson Solver on Octree in Parallel Computing
The library currently developed by the Ewha University team includes algorithms for Poisson calculations on arbitrary 2D quadtree and 3D octree structures, Hodge-Helmholtz decomposition for fluid, and fluid-solid interaction (FSI) calculations for incompressible fluids and rigid bodies, implemented in C language. Many natural and social phenomena are neither uniform nor consistent in their spatial and temporal scales. This creates a significant demand for quadtree/octree libraries that facilitate efficient scale distribution while being easy and fast to implement. We aim to deploy this library on a supercomputer to solve large-scale FSI problems as well as fluid, electromagnetics, and architectural problems related to Poisson equations. The final implementation of this library is planned to use the Modified Sparse Row (MSR) format, the standard for storing sparse or semi-sparse matrices. This allows us to leverage existing LAPACK/BLAS or comparable linear solvers embedded in supercomputers to efficiently handle linear problems.
## Poisson Solver

The traditional Poisson solver on OcTree was developed for research purposes to verify superconvergence and stability, yet it had structural inefficiencies and slow computation speed. In this module, the OcTree structure has been efficiently restructured to implement the second-order convergence code of Losasso, based on the latest research from Stanford University. This enables rapid OcTree generation in general domains, overcoming the challenging T-junction issue in OcTree by implementing Losasso’s method. Additionally, recent research by our team was the first to prove that solutions obtained using Losasso’s method achieve second-order convergence, with derivatives converging at 1.5-order.

## Incompressible Fluid

The developed Poisson solver on OcTree can be utilized as a numerical solver for various equations, particularly for calculating fluid pressure in the projection method used for incompressible fluids. The developed Poisson solver on OcTree, based on a MAC grid, naturally provides symmetry between the gradient and divergence. This symmetry enabled us to achieve L2 stability of the Hodge decomposition in the projection step of the incompressible fluid solver.



## Fluid-Solid Interaction

The FSI module addresses the movement analysis of a system where incompressible fluid interacts with rigid bodies. Due to its complexity and multiscale properties, introducing a multiscale grid like OcTree is crucial for tackling this problem. Traditional methods for computing FSI have used a staggered approach, where the fluid and solid movements are computed sequentially while the other remains fixed.

In this module, however, a monolithic method is employed, which simultaneously solves fluid and solid interactions rather than in a staggered manner. Although more complex, the monolithic method has the advantage of preserving kinetic energy, and the rigorous mathematical analysis proposed by our team clearly guarantees a convergence rate of 1.5 for accuracy.
