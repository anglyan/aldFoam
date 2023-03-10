/*---------------------------------------------------------------------------

Copyright 2019 Argonne UChicago LLC

This file is part of aldFoam.

    aldFoam is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    aldFoam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with aldFoam.  If not, see <https://www.gnu.org/licenses/>.


aldFoam solves a time-dependent reactive transport of a precursor inside
a reactor subject to two types of heterogeneous processes: a self-limited
pathway characteristic of atomic layer deposition and a non-self limited
loss that can correspond to either a CVD or recombination process.

Application
    aldFoam

Description
    Solves the reactive transport of a gaseous species under self-limited and
    non self-limited surface reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createSurfaceChem.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix N1Eqn
            (
                fvm::ddt(N1)
              + fvm::div(phi, N1)
              - fvm::laplacian(D1, N1)
             ==
                fvModels.source(N1)
            );

            N1Eqn.relax();
            fvConstraints.constrain(N1Eqn);
            N1Eqn.solve();
            fvConstraints.constrain(N1);

        }


        #include "solveSurfaceChem.H"

        runTime.write();


    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
