# What is aldFoam

aldFoam is a package to solve the reactive transport under self-limited
processes such as Atomic Layer Deposition (ALD) and Atomic Layer Etching (ALE).
aldFoam comprises a number of solvers and boundary conditions.
aldFoam is built on top of [OpenFOAM](https://openfoam.org/), a C++ toolbox
to solve partial differential equations using finite volume methods.

The version in the master branch is meant to work with OpenFoam-10.

# Install instructions

In order to install aldFoam you need to have OpenFOAM installed in your
system. From github:

```shell
git clone https://github.com/aldsim/aldFoam.git
cd aldFoam
cd firstorder
wmake
cd ..
cd dose
wmake
cd ..
cd aldFoam
wmake
```

This gives you a basic installer to simulate the reactive transport
of a single ALD or ALE precursor.

# Publication

A description of aldFoam can be found in the manuscript:

A. Yanguas-Gil, J. A. Libera, and J. W. Elam, Reactor scale simulations of ALD and ALE: ideal and non-ideal self-limited processes in a cylindrical and a 300 mm wafer cross-flow reactor, [arXiv:2106.07132](https://arxiv.org/abs/2106.07132)

The manuscript was published in the 
[Journal of Vacuum Science & Technology A](https://doi.org/10.1116/6.0001212)

# Authors

aldFoam was developed at Argonne National Laboratory. The following
scientist have been involved in the project:

  * Angel Yanguas-Gil, <ayg@anl.gov>, Lead and founder.
  * Jeffrey W Elam

# Release info

aldFoam currently contains the core solver presented in the manuscript. Other
effects, including the presence of competing byproducts and soft-saturating processes will be added as separate solvers in the near future.

Due to changes in OpenFOAM aldFoam does not work with newer versions. Due to
OpenFOAM lack of backwards compatibility in the creation of custom solvers,
it is unclear if we will be able to support newer versions.

# Copyright and license

Copyright (2019-2025) UChicago Argonne, LLC

aldFoam is distributed under the terms of the GPLv3 license. A copy
of the license is included in this distribution.
