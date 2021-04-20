/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kineticEnergy.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(kineticEnergy, 0);
    addToRunTimeSelectionTable(functionObject, kineticEnergy, dictionary);
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::kineticEnergy::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "kineticEnergy ");
    writeCommented(os, "Time");
    writeTabbed(os, "fieldIntegrate");
    os << endl;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::kineticEnergy::calc()
{
    if
        (
            foundObject<fluidThermo>(fluidThermo::dictName)
            )
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = lookupObject<volVectorField>("U");

        volScalarField V
        (
            IOobject
            (
                mesh_.V().name(),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(mesh_.V().dimensions(), Zero),
            calculatedFvPatchField<scalar>::typeName
        );

        V.ref() = mesh_.V();

        return sum(0.5 * mag(U) * mag(U) * thermo.rho() * V());
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::kineticEnergy::kineticEnergy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict)
{
    read(dict);
    writeFileHeader(file());

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::kineticEnergy::~kineticEnergy()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::functionObjects::kineticEnergy::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    //field_ = word(dict.lookup("field"));


    return true;
}


bool Foam::functionObjects::kineticEnergy::execute()
{
    return true;
}


bool Foam::functionObjects::kineticEnergy::end()
{
    return true;
}


bool Foam::functionObjects::kineticEnergy::write()
{
    const scalar intValue = calc();


    //Log << type() << " " << kineticEnergy.name() << " is " << intValue << endl;
    Log << type() << " " << "Kinetic Energy is " << intValue << endl;


    if (Pstream::master())
    {
        writeTime(file());

        file()
            << token::TAB << intValue
            << endl;
    }



    return true;
}


// ************************************************************************* //