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



Info<< "Reading field N1\n" << endl;

volScalarField N1
(
    IOobject
    (
        "N1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field Exp1\n" << endl;

volScalarField Exp1
(
    IOobject
    (
        "Exp1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field cov\n" << endl;

volScalarField cov
(
    IOobject
    (
        "cov",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field beta\n" << endl;

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field growth\n" << endl;

volScalarField growth
(
    IOobject
    (
        "growth",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Nold = N1;

Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity D1\n" << endl;

dimensionedScalar D1
(
    transportProperties.lookup("D1")
);

dimensionedScalar vth
(
    transportProperties.lookup("vth")
);

#include "createPhi.H"

#include "createFvModels.H"
#include "createFvConstraints.H"
