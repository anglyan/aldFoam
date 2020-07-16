# What is aldFoam

aldFoam is a package to solve the reactive transport under self-limited
processes such as Atomic Layer Deposition (ALD) and Atomic Layer Etching (ALE).
aldFoam comprises a number of solvers and boundary conditions.
aldFoam is built on top of [OpenFOAM](https://openfoam.org/), a C++ toolbox
to solve partial differential equations using finite volume methods.

# Disclaimer

What you are seeing here is a pre-released version of aldFoam. The
first public release will happen in late July 2020.

# Install instructions

In order to install aldFoam you need to have OpenFOAM installed in your
system. From github:

```shell
git clone https://github.com/anglyan/aldFoam.git
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

# Authors

aldFoam was developed at Argonne National Laboratory. The following
scientist have been involved in the project:

  * Angel Yanguas-Gil, <ayg@anl.gov>, Lead and founder.
  * Jeffrey W Elam

# Copyright and license

Copyright (2019) UChicago Argonne, LLC

aldFoam is distributed under the terms of the GPLv3 license. A copy
of the license is included in this distribution.
