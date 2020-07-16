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

Info << "Reading process data from processData\n" << endl;

IOdictionary processData
(
    IOobject
    (
        "processData",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar pressure = readScalar(processData.lookup("pressure"));
scalar temperature = readScalar(processData.lookup("temperature"));


Info << "Reading surface processes.\n" << endl;

IOdictionary surfaceProcesses

(
    IOobject
    (
        "surfaceProcesses",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info << "Reading surface chemistry patches from surfacePatches" << endl;

label npatches = readLabel(surfaceProcesses.lookup("numPatches"));

dictionary patchKeys(surfaceProcesses.subDict("surfacePatches"));
wordList surfacePatches(npatches);
scalarList beta0(npatches);
scalarList beta1(npatches);
scalarList betarec(npatches);
scalarList nu0(npatches);


forAll(surfacePatches, i){
    surfacePatches[i] = (word) patchKeys.keys()[i];
    dictionary surfchem(patchKeys.subDict(surfacePatches[i]));
    scalar sitearea = readScalar(surfchem.lookup("sitearea"));
    nu0[i] = 0.25*vth.value()*pressure/(1.38e-23*temperature)*sitearea;
    beta0[i] = readScalar(surfchem.lookup("beta0"));
    beta1[i] = readScalar(surfchem.lookup("beta1"));
    betarec[i] = readScalar(surfchem.lookup("betarec"));
    label patchID =
        mesh.boundaryMesh().findPatchID(surfacePatches[i]);

    scalarField & betap = beta.boundaryFieldRef()[patchID];
    betap = beta0[i] + beta1[i] + betarec[i];
    Info << surfacePatches[i] << endl;
    Info << nu0[i] << endl;
}
