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


\*---------------------------------------------------------------------------*/


        dimensionedScalar dt = runTime.deltaT();

        forAll(surfacePatches, i){
            label patchID =
                mesh.boundaryMesh().findPatchID(surfacePatches[i]);
            Info << "Patch ID:" << patchID << endl;

            scalarField & N1b = N1.boundaryFieldRef()[patchID];
            scalarField & Expb = Exp1.boundaryFieldRef()[patchID];
            scalarField & Noldb = Nold.boundaryFieldRef()[patchID];
            scalarField & growthp = growth.boundaryFieldRef()[patchID];


            Expb =  nu0[i]*N1b*dt.value();

            scalarField & betap = beta.boundaryFieldRef()[patchID];
            scalarField & covp = cov.boundaryFieldRef()[patchID];
            covp = (covp + beta0[i]*Expb)/(1+beta0[i]*Expb);
            betap = beta0[i]*(1-covp)+beta1[i]+betarec[i];
            growthp += betap*Expb;
        }
